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
using base::not_null;
using geometry::Instant;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::RotatingBody;
using quantities::Time;

class OrbitAnalyser {
 public:
  struct Parameters {
    Ephemeris<Barycentric>::Guard guard;
    Instant first_time;
    DegreesOfFreedom<Barycentric> first_degrees_of_freedom;
    Time mission_duration;
    not_null<RotatingBody<Barycentric> const*> primary;
    bool shutdown = false;
  };

  struct Analysis {
    Instant first_time;
    Time mission_duration;
    not_null<RotatingBody<Barycentric> const*> primary;
    std::optional<OrbitalElements> elements;
    std::optional<OrbitGroundTrack> ground_track;
  };

  OrbitAnalyser(
      not_null<Ephemeris<Barycentric>*> ephemeris,
      Ephemeris<Barycentric>::FixedStepParameters analysed_trajectory_parameters);

 private:
  void RepeatedlyAnalyseOrbit();

  not_null<Ephemeris<Barycentric>*> const ephemeris_;
  Ephemeris<Barycentric>::FixedStepParameters const
      analysed_trajectory_parameters_;

  mutable absl::Mutex lock_;
  std::thread analyser_;
  std::optional<Analysis> analysis_;
  std::optional<Parameters> parameters_ GUARDED_BY(lock_);
  std::optional<Analysis> next_analysis_ GUARDED_BY(lock_);
  std::atomic_int8_t next_analysis_percentage_;
};

}  // namespace internal_orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
