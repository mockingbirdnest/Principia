#pragma once

#include <atomic>
#include <optional>
#include <thread>

#include "absl/synchronization/mutex.h"
#include "astronomy/orbit_ground_track.hpp"
#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "base/jthread.hpp"
#include "base/not_null.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _orbit_analyser {
namespace internal {

using namespace principia::astronomy::_orbit_ground_track;
using namespace principia::astronomy::_orbit_recurrence;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::base::_jthread;
using namespace principia::base::_not_null;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_named_quantities;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_rotating_body;
using namespace principia::quantities::_quantities;

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
    RotatingBody<Barycentric> const* primary() const;
    std::optional<Interval<Length>> radial_distance_interval() const;
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
    explicit Analysis(Instant const& first_time);

    Instant first_time_;
    Time mission_duration_;
    RotatingBody<Barycentric> const* primary_ = nullptr;
    std::optional<Interval<Length>> radial_distance_interval_;
    std::optional<OrbitalElements> elements_;
    std::optional<OrbitRecurrence> closest_recurrence_;
    std::optional<OrbitRecurrence> recurrence_;
    std::optional<OrbitGroundTrack> ground_track_;
    std::optional<OrbitGroundTrack::EquatorCrossingLongitudes>
        equatorial_crossings_;

    friend class OrbitAnalyser;
  };

  struct Parameters {
    Instant first_time;
    DegreesOfFreedom<Barycentric> first_degrees_of_freedom;
    Time mission_duration;
    // The analyser may compute the trajectory up to |extended_mission_duration|
    // to ensure that at least one revolution is analysed.
    std::optional<Time> extended_mission_duration;
  };

  OrbitAnalyser(not_null<Ephemeris<Barycentric>*> ephemeris,
                Ephemeris<Barycentric>::FixedStepParameters
                    analysed_trajectory_parameters);

  virtual ~OrbitAnalyser();

  // Cancel any computation in progress, causing the next call to
  // |RequestAnalysis| to be processed as fast as possible.
  void Interrupt();

  // Sets the parameters that will be used for the computation of the next
  // analysis.
  void RequestAnalysis(Parameters const& parameters);

  // The last value passed to |RequestAnalysis|.
  std::optional<Parameters> const& last_parameters() const;

  // Sets |analysis()| to the latest computed analysis.
  void RefreshAnalysis();

  // Mutable so that the caller can call |SetRecurrence| and |ResetRecurrence|.
  Analysis* analysis();

  // The result is in [0, 1]; it tracks the progress of the computation of the
  // next analysis.  Note that a new analysis may be ready even if this is not
  // equal to 1, if the analyser is working on a subsequent request.
  double progress_of_next_analysis() const;

 private:
  absl::Status AnalyseOrbit(Parameters const& parameters);

  // Flows the |trajectory| with a fixed step integrator using the given
  // |parameters|.  This is done in small increments and
  // |progress_of_next_analysis_| is updated after each increment to be able to
  // display a progress bar.  Returns the actual mission duration.
  absl::StatusOr<Time> FlowWithProgressBar(
      Parameters const& parameters,
      Time const& smallest_osculating_period,
      DiscreteTrajectory<Barycentric>& trajectory);

  not_null<Ephemeris<Barycentric>*> const ephemeris_;
  Ephemeris<Barycentric>::FixedStepParameters const
      analysed_trajectory_parameters_;

  std::optional<Parameters> last_parameters_;

  std::optional<Analysis> analysis_;

  mutable absl::Mutex lock_;
  jthread analyser_;
  // The |analyser_| is idle:
  // — if it is not joinable, e.g. because it was stopped by |Interrupt()|, or
  // — if it is done computing |next_analysis_| and has stopped or is about to
  //   stop executing.
  // If it is joined once idle (and joinable), it will not attempt to acquire
  // |lock_|.
  bool analyser_idle_ GUARDED_BY(lock_) = true;
  // |next_analysis_| is set by the |analyser_| thread; it is read and cleared
  // by the main thread.
  std::optional<Analysis> next_analysis_ GUARDED_BY(lock_);
  // |progress_of_next_analysis_| is set by the |analyser_| thread; it tracks
  // progress in computing |next_analysis_|.
  std::atomic<double> progress_of_next_analysis_ = 0;
};

}  // namespace internal

using internal::OrbitAnalyser;

}  // namespace _orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia

namespace principia::ksp_plugin {
using namespace principia::ksp_plugin::_orbit_analyser;
}  // namespace principia::ksp_plugin
