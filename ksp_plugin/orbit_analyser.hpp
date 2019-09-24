#pragma once

#include <atomic>
#include <optional>
#include <thread>

#include "absl/synchronization/mutex.h"
#include "astronomy/orbit_ground_track.hpp"
#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_orbit_analyser {

using astronomy::OrbitalElements;
using astronomy::OrbitGroundTrack;
using astronomy::OrbitRecurrence;
using base::not_null;
using geometry::Instant;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::RotatingBody;
using quantities::Time;

// The |OrbitAnalyser| asynchronously integrates a trajectory, and computes
// orbital elements, recurrence, and ground track properties of the resulting
// orbit.
class OrbitAnalyser {
 public:
  // The analysis stores the computed orbital characteristics.  It is publicly
  // mutable via |SetRecurrence| and |ResetRecurrence| to allow the caller to
  // consider a nominal recurrence other than the one deduced from the orbital
  // elements: analysing the precomputed ground track with respect to a
  // different recurrence is relatively cheap, so it is inconvenient to wait for
  // a whole new analysis to do so, but doing it at every frame is still
  // wasteful, so we cache that in the |Analysis|.
  class Analysis {
   public:
    Instant const& first_time() const;
    Time const& mission_duration() const;
    RotatingBody<Barycentric> const& primary() const;
    std::optional<OrbitalElements> const& elements() const;
    std::optional<OrbitRecurrence> const& recurrence() const;
    std::optional<OrbitGroundTrack> const& ground_track() const;
    // |equatorial_crossings().has_value()| if and only if
    // |recurrence().has_value && ground_track().has_value()|;
    // |*equatorial_crossings()| is
    //   ground_track()->equator_crossing_longitudes(
    //       *recurrence(), /*first_ascending_pass_index=*/1)
    // precomputed to avoid performing this calculation at every frame.
    std::optional<OrbitGroundTrack::EquatorCrossingLongitudes> const&
    equatorial_crossings() const;

    // Sets |recurrence|, updating |equatorial_crossings| if needed.
    void SetRecurrence(OrbitRecurrence const& recurrence);
    // Resets |recurrence| to a value deduced from |*elements| by
    // |OrbitRecurrence::ClosestRecurrence|, or to nullopt if
    // |!elements.has_value()|, updating |equatorial_crossings| if needed.
    void ResetRecurrence();

   private:
    Analysis(Instant const& first_time,
             not_null<RotatingBody<Barycentric> const*> primary);

    Instant first_time_;
    Time mission_duration_;
    not_null<RotatingBody<Barycentric> const*> primary_;
    std::optional<OrbitalElements> elements_;
    std::optional<OrbitRecurrence> closest_recurrence_;
    std::optional<OrbitRecurrence> recurrence_;
    std::optional<OrbitGroundTrack> ground_track_;
    std::optional<OrbitGroundTrack::EquatorCrossingLongitudes>
        equatorial_crossings_;

    friend class OrbitAnalyser;
  };

  OrbitAnalyser(not_null<Ephemeris<Barycentric>*> ephemeris,
                Ephemeris<Barycentric>::FixedStepParameters const&
                    analysed_trajectory_parameters);
  ~OrbitAnalyser();

  // Sets the parameters that will be used for the computation of the next analysis.
  void RequestAnalysis(
      Instant const& first_time,
      DegreesOfFreedom<Barycentric> const& first_degrees_of_freedom,
      Time const& mission_duration,
      not_null<RotatingBody<Barycentric> const*> primary);

  // Sets |analysis()| to the latest computed analysis.
  void RefreshAnalysis();

  // Mutable so that the caller can call |SetRecurrence| and |ResetRecurrence|.
  Analysis* analysis();

  // The result is in [0, 1]; it tracks the progress of the computation of the
  // next analysis.  Note that a new analysis may be ready even if this is not
  // equal to 1, if a subsequent request is being processed.
  double progress_of_next_analysis() const;

 private:
  struct Parameters {
    Ephemeris<Barycentric>::Guard guard;
    Instant first_time;
    DegreesOfFreedom<Barycentric> first_degrees_of_freedom;
    Time mission_duration;
    not_null<RotatingBody<Barycentric> const*> primary;
  };

  void RepeatedlyAnalyseOrbit();

  not_null<Ephemeris<Barycentric>*> const ephemeris_;
  Ephemeris<Barycentric>::FixedStepParameters const
      analysed_trajectory_parameters_;

  std::optional<Analysis> analysis_;

  mutable absl::Mutex lock_;
  std::thread analyser_;
  // |parameters_| is set by the main thread; it is read and cleared by the
  // |analyser_| thread.
  std::optional<Parameters> parameters_ GUARDED_BY(lock_);
  // |next_analysis_| is set by the |analyser_| thread; it is read and cleared
  // by the main thread.
  std::optional<Analysis> next_analysis_ GUARDED_BY(lock_);
  // |progress_of_next_analysis_| is set by the |analyser_| thread; it tracks
  // progress in computing |next_analysis_|.
  std::atomic<double> progress_of_next_analysis_ = 0;
  // |keep_analysing_| is tested by the |analyser_| thread, which cooperatively
  // aborts if it is false; it is set at construction, and cleared by the main
  // thread at destruction.
  std::atomic_bool keep_analysing_ = true;
};

}  // namespace internal_orbit_analyser

using internal_orbit_analyser::OrbitAnalyser;

}  // namespace ksp_plugin
}  // namespace principia
