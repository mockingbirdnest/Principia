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

class OrbitAnalyser {
 public:

  struct Analysis {
    Instant first_time;
    Time mission_duration;
    not_null<RotatingBody<Barycentric> const*> primary;
    std::optional<OrbitalElements> elements;
    std::optional<OrbitRecurrence> auto_detected_recurrence;
    std::optional<OrbitRecurrence> chosen_recurrence;
    std::optional<OrbitGroundTrack> ground_track;
    std::optional<OrbitGroundTrack::EquatorCrossingLongitudes>
        equator_crossings_for_chosen_recurrence;
  };

  OrbitAnalyser(
      not_null<Ephemeris<Barycentric>*> ephemeris,
      Ephemeris<Barycentric>::FixedStepParameters analysed_trajectory_parameters);
  ~OrbitAnalyser();

  // Sets the parameters that will be used for the computation of the next analysis.
  void RequestAnalysis(
      Instant const& first_time,
      DegreesOfFreedom<Barycentric> const& first_degrees_of_freedom,
      Time const& mission_duration,
      not_null<RotatingBody<Barycentric> const*> primary);

  // Sets |analysis()| to the latest computed analysis.
  void RefreshAnalysis();

  std::optional<Analysis>& analysis();

  int8_t next_analysis_percentage() const;

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
  // |next_analysis_percentage_| is set by the |analyser_| thread; it tracks
  // progress in computing |next_anlaysis_|.
  std::atomic_int8_t next_analysis_percentage_ = 0;
  // |keep_analysing_| is tested |analyser_| thread, which cooperatively aborts
  // if it is false; it is set at construction, and cleared by the main thread
  // at destruction.
  std::atomic_bool keep_analysing_ = true;
};

}  // namespace internal_orbit_analyser

using internal_orbit_analyser::OrbitAnalyser;

}  // namespace ksp_plugin
}  // namespace principia
