#pragma once

#include <optional>

#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "geometry/interval.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/rotating_body.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_ground_track {

using geometry::Instant;
using geometry::Interval;
using physics::DiscreteTrajectory;
using physics::RotatingBody;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Time;

class OrbitGroundTrack {
 public:
  struct MeanSun {
    Instant epoch;
    // This mean longitude of the sun is with respect the axes of
    // |PrimaryCentred|.
    Angle mean_longitude_at_epoch;
    // This year is the period of the above mean longitude; which kind of year
    // it is depends on the definition of |PrimaryCentred|.  If |PrimaryCentred|
    // has ICRS axes, it is the sidereal year; if |PrimaryCentred| is the
    // reference frame of the equator & equinox of the date, it is the tropical
    // year.
    Time year;
  };

  template<typename PrimaryCentred, typename Inertial>
  static OrbitGroundTrack ForTrajectory(
      DiscreteTrajectory<PrimaryCentred> const& trajectory,
      RotatingBody<Inertial> const& primary,
      std::optional<OrbitRecurrence> const& nominal_recurrence,
      std::optional<MeanSun> const& mean_sun);

  // The interval spanned by the longitudes of the crossings of the equator on
  // the ascending passes, reduced in the sense that they are offset to
  // compensate for the nominal equatorial shift, with the initial value reduced
  // to an eastward grid interval from the zero meridian (longitudes [0, δ]).
  // This is populated only if a nominal recurrence was provided.
  std::optional<Interval<Angle>> const&
  reduced_longitudes_of_equator_crossings_of_ascending_passes() const;

  // The interval spanned by the local mean solar times at the ascending nodes.
  // This is populated only if a mean sun was provided.
  // The initial value lies in [0, 2π], with π being noon.
  std::optional<Interval<Angle>> const& mean_solar_times_of_ascending_nodes()
      const;
  // Same as above with the descending nodes.
  std::optional<Interval<Angle>> const& mean_solar_times_of_descending_nodes()
      const;

 private:
  OrbitGroundTrack() = default;

  std::optional<Interval<Angle>>
      reduced_longitudes_of_equator_crossings_of_ascending_passes_;
  std::optional<Interval<Angle>> mean_solar_times_of_ascending_nodes_;
  std::optional<Interval<Angle>> mean_solar_times_of_descending_nodes_;
};

}  // namespace internal_orbit_ground_track

using internal_orbit_ground_track::OrbitGroundTrack;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_ground_track_body.hpp"
