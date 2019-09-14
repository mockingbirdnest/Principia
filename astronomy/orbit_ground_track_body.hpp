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

// The resulting angles are neither normalized nor unwound.
template<typename PrimaryCentred, typename Inertial>
std::vector<Angle> PlanetocentricLongitudes(
    DiscreteTrajectory<PrimaryCentred> const& nodes,
    RotatingBody<Inertial> const& primary) {
  std::vector<Angle> longitudes;
  for (auto const& node : nodes) {
    longitudes.push_back(
        CelestialLongitude(node.degrees_of_freedom.position()) -
        primary.AngleAt(node.time) - π / 2 * Radian);
  }
  return longitudes;
}

// The resulting angle is not normalized.
template<typename Iterator>
Angle MeanSolarTime(Iterator const& it,
                    OrbitGroundTrack::MeanSun const& mean_sun) {
  Time const t = it->time - mean_sun.epoch;
  return π * Radian + CelestialLongitude(it->degrees_of_freedom.position()) -
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

inline Interval<Angle> ReducedLongitudesOfEquatorialCrossings(
    std::vector<Angle> const& longitudes_of_equatorial_crossings,
    OrbitRecurrence const& nominal_recurrence,
    Angle const& initial_offset) {
  Interval<Angle> reduced_longitudes;
  std::optional<Angle> reduced_longitude;
  for (int n = 0; n != longitudes_of_equatorial_crossings.size(); ++n) {
    Angle const longitude = longitudes_of_equatorial_crossings[n];
    Angle const offset_longitude =
        longitude - initial_offset - n * nominal_recurrence.equatorial_shift();
    reduced_longitude = reduced_longitude.has_value()
                            ? UnwindFrom(*reduced_longitude, offset_longitude)
                            : Mod(offset_longitude, 2 * π * Radian);
    reduced_longitudes.Include(*reduced_longitude);
  }
  return reduced_longitudes;
}

Interval<Angle>
OrbitGroundTrack::EquatorCrossingLongitudes::longitudes_reduced_to_pass(
    int const pass_index) const {
  if (pass_index % 2 == 1) {
    Angle const shift = ((pass_index - 1) / 2) *
                        nominal_recurrence_.equatorial_shift();
    Angle const min = Mod(ascending_longitudes_reduced_to_pass_1_.min + shift,
                          2 * π * Radian);
    Angle const max =
        UnwindFrom(min, ascending_longitudes_reduced_to_pass_1_.max + shift);
    return {min, max};
  } else {
    Angle const shift =
        ((pass_index - 2) / 2) * nominal_recurrence_.equatorial_shift();
    Angle const min = Mod(descending_longitudes_reduced_to_pass_2_.min + shift,
                          2 * π * Radian);
    Angle const max =
        UnwindFrom(min, descending_longitudes_reduced_to_pass_2_.max + shift);
    return {min, max};
  }
}

inline OrbitGroundTrack::EquatorCrossingLongitudes::EquatorCrossingLongitudes(
    OrbitRecurrence const& nominal_recurrence)
    : nominal_recurrence_(nominal_recurrence) {}

template<typename PrimaryCentred, typename Inertial>
OrbitGroundTrack OrbitGroundTrack::ForTrajectory(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    RotatingBody<Inertial> const& primary,
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
  ground_track.longitudes_of_equator_crossings_of_ascending_passes_ =
      PlanetocentricLongitudes(ascending_nodes, primary);
  ground_track.longitudes_of_equator_crossings_of_descending_passes_ =
      PlanetocentricLongitudes(descending_nodes, primary);
  return ground_track;
}

inline OrbitGroundTrack::EquatorCrossingLongitudes
OrbitGroundTrack::equator_crossing_longitudes(
    OrbitRecurrence const& nominal_recurrence,
    int const first_ascending_pass_index) const {
  EquatorCrossingLongitudes equator_crossings{nominal_recurrence};
  if (longitudes_of_equator_crossings_of_ascending_passes_.empty()) {
    return equator_crossings;
  }

  Angle const initial_offset =
      Mod(longitudes_of_equator_crossings_of_ascending_passes_.front(),
          nominal_recurrence.grid_interval()) -
      ((first_ascending_pass_index - 1) / 2) *
          nominal_recurrence.equatorial_shift();
  equator_crossings.ascending_longitudes_reduced_to_pass_1_ =
      ReducedLongitudesOfEquatorialCrossings(
          longitudes_of_equator_crossings_of_ascending_passes_,
          nominal_recurrence,
          initial_offset);
  equator_crossings.descending_longitudes_reduced_to_pass_2_ =
      ReducedLongitudesOfEquatorialCrossings(
          longitudes_of_equator_crossings_of_descending_passes_,
          nominal_recurrence,
          initial_offset);
  return equator_crossings;
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
