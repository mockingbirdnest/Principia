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

class OrbitGroundTrack {
 public:
  template<typename PrimaryCentred, typename Inertial>
  OrbitGroundTrack ForTrajectory(
      DiscreteTrajectory<PrimaryCentred> const& trajectory,
      RotatingBody<Inertial> const& primary,
      std::optional<OrbitRecurrence> const& nominal_recurrence);

  // The interval spanned by the geographical longitudes of the ascending nodes,
  // compensating for the nominal equatorial shift, with the initial value
  // reduced to an eastward grid interval (longitudes [0, δ[).
  // This is populated only if a recurrence was provided.
  std::optional<Interval<Angle>> const& reduced_longitude_of_ascending_node()
      const;

  // The intervals spanned by the geographical coordinates of the apoapsides.
  // These are populated only if a recurrence was provided.
  // The longitude is reduced in the same way as
  // `reduced_longitude_of_ascending_node`.
  std::optional<Interval<Angle>> const& reduced_longitude_of_apoapsis() const;
  std::optional<Interval<Angle>> const& latitude_of_apoapsis() const;

  // The interval spanned by the local mean solar times at the ascending nodes.
  // This is populated only if a mean sun was provided.
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
