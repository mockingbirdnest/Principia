#pragma once

#include <optional>

#include "astronomy/orbit_recurrence.hpp"
#include "geometry/interval.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/rotating_body.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_ground_track {

using geometry::Interval;
using physics::DiscreteTrajectory;
using physics::RotatingBody;
using quantities::Angle;

class OrbitGroundTrack {
 public:
  template<typename PrimaryCentred, typename Inertial>
  static OrbitGroundTrack ForTrajectory(
      DiscreteTrajectory<PrimaryCentred> const& trajectory,
      RotatingBody<Inertial> const& primary,
      std::optional<OrbitRecurrence> const& nominal_recurrence);

  // The interval spanned by the geographical longitudes of the ascending nodes,
  // compensating for the nominal equatorial shift, with the initial value
  // reduced to an eastward grid interval (longitudes [0, δ]).
  // This is populated only if a nominal recurrence was provided.
  std::optional<Interval<Angle>> const& reduced_longitude_of_ascending_node()
      const;

  // The time average of each geographical coordinate of the ground track.
  // These are populated only if the nominal recurrence is [1, 0, 1]
  // (synchronous orbit).
  std::optional<Angle> const& average_longitude() const;
  std::optional<Angle> const& average_latitude() const;

  // The interval spanned by the local mean solar times at the ascending nodes.
  // This is populated only if a mean sun was provided.
  // The initial value lies in [-π, π], with 0 being noon.
  // TODO(egg): since we do not have a mean sun class at this time, this is
  // always nullopt.
  std::optional<Interval<Angle>> const& mean_solar_time_of_ascending_node()
      const;

 private:
  OrbitGroundTrack() = default;

  std::optional<Interval<Angle>> reduced_longitude_of_ascending_node_;
  std::optional<Interval<Angle>> reduced_longitude_of_apoapsis_;
  std::optional<Interval<Angle>> latitude_of_apoapsis_;
  std::optional<Interval<Angle>> mean_solar_time_of_ascending_node_;
};

}  // namespace internal_orbit_ground_track

using internal_orbit_ground_track::OrbitGroundTrack;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_ground_track_body.hpp"
