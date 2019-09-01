#pragma once

#include "astronomy/orbit_ground_track.hpp"

#include <limits>

#include "physics/apsides.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_ground_track {

using geometry::Position;
using geometry::Vector;
using physics::ComputeNodes;
using quantities::Mod;
using quantities::UnwindFrom;
using quantities::si::Radian;

// Note that the origin of this celestial longitude is arbitrary: it is not the
// node of the orbit around the sun (the equinox).  If |PrimaryCentred| is the
// GCRS, this is the right ascension (with respect to the mean equinox of
// J2000).
template<typename PrimaryCentred>
Angle CelestialLongitude(Position<PrimaryCentred> const& q) {
  // TODO(egg): |ToSpherical| is wasteful, as we discard the latitude.
  return (q - PrimaryCentred::origin).coordinates().ToSpherical().longitude;
}

// The resulting angle is not normalized.
template<typename Iterator, typename Inertial>
Angle PlanetocentricLongitude(Iterator const& it,
                              RotatingBody<Inertial> const& primary) {
  return CelestialLongitude(it.degrees_of_freedom().position()) -
         primary.AngleAt(it.time()) - π / 2 * Radian;
}

// The resulting angle is not normalized.
template<typename Iterator>
Angle MeanSolarTime(Iterator const& it,
                    OrbitGroundTrack::MeanSun const& mean_sun) {
  Time const t = it.time() - mean_sun.epoch;
  return π * Radian + CelestialLongitude(it.degrees_of_freedom().position()) -
         (mean_sun.mean_longitude_at_epoch +
          (2 * π * Radian * t / mean_sun.year));
}

template<typename PrimaryCentred>
Interval<Angle> MeanSolarTimesOfNodes(
    DiscreteTrajectory<PrimaryCentred> const& nodes,
    OrbitGroundTrack::MeanSun const& mean_sun) {
  Interval<Angle> mean_solar_times;
  std::optional<Angle> mean_solar_time;
  for (auto node = nodes.begin(); node != nodes.end(); ++node) {
    if (mean_solar_time.has_value()) {
      mean_solar_time =
          UnwindFrom(*mean_solar_time, MeanSolarTime(node, mean_sun));
    } else {
      mean_solar_time = Mod(MeanSolarTime(node, mean_sun), 2 * π * Radian);
    }
    mean_solar_times.Include(*mean_solar_time);
  }
  return mean_solar_times;
}

template<typename PrimaryCentred, typename Inertial>
Interval<Angle> ReducedLongitudesOfEquatorialCrossings(
    DiscreteTrajectory<PrimaryCentred> const& nodes,
    RotatingBody<Inertial> const& primary,
    OrbitRecurrence const& nominal_recurrence,
    std::optional<Angle>& initial_offset) {
  Interval<Angle> reduced_longitudes;
  std::optional<Angle> reduced_longitude;
  for (auto [node, n] = std::make_pair(nodes.begin(), 0);  // NOLINT
       node != nodes.end();
       ++node, ++n) {
    Angle const planetocentric_longitude =
        PlanetocentricLongitude(node, primary);

    if (!initial_offset.has_value()) {
      reduced_longitude =
          Mod(planetocentric_longitude, nominal_recurrence.grid_interval());
      initial_offset = planetocentric_longitude - *reduced_longitude;
    } else {
      Angle const offset_longitude = planetocentric_longitude -
                                     *initial_offset -
                                     n * nominal_recurrence.equatorial_shift();
      reduced_longitude = reduced_longitude.has_value()
                              ? UnwindFrom(*reduced_longitude, offset_longitude)
                              : Mod(offset_longitude, 2 * π * Radian);
    }
    reduced_longitudes.Include(*reduced_longitude);
  }
  return reduced_longitudes;
}

template<typename PrimaryCentred, typename Inertial>
OrbitGroundTrack OrbitGroundTrack::ForTrajectory(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    RotatingBody<Inertial> const& primary,
    std::optional<OrbitRecurrence> const& nominal_recurrence,
    std::optional<MeanSun> const& mean_sun) {
  DiscreteTrajectory<PrimaryCentred> ascending_nodes;
  DiscreteTrajectory<PrimaryCentred> descending_nodes;
  OrbitGroundTrack ground_track;
  ComputeNodes(trajectory.begin(),
               trajectory.end(),
               Vector<double, PrimaryCentred>({0, 0, 1}),
               /*max_points=*/std::numeric_limits<int>::max(),
               ascending_nodes,
               descending_nodes);
  if (mean_sun.has_value()) {
    if (!ascending_nodes.Empty()) {
      ground_track.mean_solar_times_of_ascending_nodes_ =
          MeanSolarTimesOfNodes(ascending_nodes, *mean_sun);
    }
    if (!descending_nodes.Empty()) {
      ground_track.mean_solar_times_of_descending_nodes_ =
          MeanSolarTimesOfNodes(descending_nodes, *mean_sun);
    }
  }
  if (nominal_recurrence.has_value()) {
    std::optional<Angle> initial_offset;
    if (!ascending_nodes.Empty()) {
      ground_track
          .reduced_longitudes_of_equator_crossings_of_ascending_passes_ =
          ReducedLongitudesOfEquatorialCrossings(
              ascending_nodes, primary, *nominal_recurrence, initial_offset);
    }
    if (!descending_nodes.Empty()) {
      ground_track
          .reduced_longitudes_of_equator_crossings_of_descending_passes_ =
          ReducedLongitudesOfEquatorialCrossings(
              descending_nodes, primary, *nominal_recurrence, initial_offset);
    }
  }
  return ground_track;
}

inline std::optional<Interval<Angle>> const&
OrbitGroundTrack::reduced_longitudes_of_equator_crossings_of_ascending_passes()
    const {
  return reduced_longitudes_of_equator_crossings_of_ascending_passes_;
}

inline std::optional<Interval<Angle>> const&
OrbitGroundTrack::reduced_longitudes_of_equator_crossings_of_descending_passes()
    const {
  return reduced_longitudes_of_equator_crossings_of_descending_passes_;
}

inline std::optional<Interval<Angle>> const&
OrbitGroundTrack::mean_solar_times_of_ascending_nodes() const {
  return mean_solar_times_of_ascending_nodes_;
}

inline std::optional<Interval<Angle>> const&
OrbitGroundTrack::mean_solar_times_of_descending_nodes() const {
  return mean_solar_times_of_descending_nodes_;
}

}  // namespace internal_orbit_ground_track
}  // namespace astronomy
}  // namespace principia
