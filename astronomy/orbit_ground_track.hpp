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
    // This mean longitude of the sun is with respect to the axes of
    // |PrimaryCentred|.
    Angle mean_longitude_at_epoch;
    // This year is the period of the above mean longitude; which kind of year
    // it is depends on the definition of |PrimaryCentred|.  If |PrimaryCentred|
    // has ICRS axes, it is the sidereal year; if |PrimaryCentred| is the
    // reference frame of the equator & equinox of the date, it is the tropical
    // year.
    Time year;
  };

  class EquatorCrossingLongitudes {
   public:
    // The longitudes of the equator crossings (ascending or descending,
    // depending on the parity of pass_index), reduced to a grid interval around
    // the pass with the given index.  This provides both an indication of how
    // well the orbit follows the nominal recurrence grid, and where that grid
    // is located.
    Interval<Angle> longitudes_reduced_to_pass(int const pass_index) const;

   private:
    EquatorCrossingLongitudes(OrbitRecurrence const& nominal_recurrence);
    OrbitRecurrence nominal_recurrence_;
    Interval<Angle> ascending_longitudes_reduced_to_pass_1_;
    Interval<Angle> descending_longitudes_reduced_to_pass_2_;
    friend class OrbitGroundTrack;
  };

  // Returns an object that describes the properties of the ground track of
  // |trajectory| as an orbit around |primary|; if |mean_sun| is provided,
  // sun-synchronicity is analysed.
  template<typename PrimaryCentred, typename Inertial>
  static OrbitGroundTrack ForTrajectory(
      DiscreteTrajectory<PrimaryCentred> const& trajectory,
      RotatingBody<Inertial> const& primary,
      std::optional<MeanSun> const& mean_sun);

  // Given a nominal recurrence and the index of the first ascending pass east
  // of the equator (which must be odd and in [1, 2 Nᴛₒ - 1], where Nᴛₒ is
  // |nominal_recurrence.number_of_revolutions()|), returns an object describing
  // how the recurrence grid is aligned in longitude, and how well the actual
  // orbit follows that grid.
  EquatorCrossingLongitudes equator_crossing_longitudes(
      OrbitRecurrence const& nominal_recurrence,
      int first_ascending_pass_index) const;

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

  std::vector<Angle> longitudes_of_equator_crossings_of_ascending_passes_;
  std::vector<Angle> longitudes_of_equator_crossings_of_descending_passes_;
  // Whether |longitudes_of_equator_crossings_of_descending_passes_.front()| is
  // the longitude of the pass preceding
  // |longitudes_of_equator_crossings_of_ascending_passes_.front()|, rather than
  // the following one.
  bool first_descending_pass_before_first_ascending_pass_;
  std::optional<Interval<Angle>> mean_solar_times_of_ascending_nodes_;
  std::optional<Interval<Angle>> mean_solar_times_of_descending_nodes_;
};

}  // namespace internal_orbit_ground_track

using internal_orbit_ground_track::OrbitGroundTrack;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_ground_track_body.hpp"
