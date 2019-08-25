#pragma once

#include "astronomy/orbit_ground_track.hpp"

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

template<typename PrimaryCentred, typename Inertial>
OrbitGroundTrack OrbitGroundTrack::ForTrajectory(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    RotatingBody<Inertial> const& primary,
    std::optional<OrbitRecurrence> const& nominal_recurrence,
    std::optional<MeanSun> const& mean_sun) {
  DiscreteTrajectory<PrimaryCentred> ascending_nodes;
  DiscreteTrajectory<PrimaryCentred> descending_nodes;
  OrbitGroundTrack ground_track;
  ComputeNodes(trajectory.Begin(),
               trajectory.End(),
               Vector<double, PrimaryCentred>({0, 0, 1}),
               /*max_points=*/std::numeric_limits<int>::max(),
               ascending_nodes,
               descending_nodes);
  if (!descending_nodes.Empty()) {
    if (mean_sun.has_value()) {
      std::optional<Angle> mean_solar_time;
      ground_track.mean_solar_times_of_descending_nodes_.emplace();
      for (auto descending_node = descending_nodes.Begin();
           descending_node != descending_nodes.End();
           ++descending_node) {
        if (mean_solar_time.has_value()) {
          mean_solar_time = UnwindFrom(
              *mean_solar_time, MeanSolarTime(descending_node, *mean_sun));
        } else {
          mean_solar_time =
              Mod(MeanSolarTime(descending_node, *mean_sun), 2 * π * Radian);
        }
        ground_track.mean_solar_times_of_descending_nodes_->Include(
            *mean_solar_time);
      }
    }
  }
  if (!ascending_nodes.Empty()) {
    if (mean_sun.has_value()) {
      std::optional<Angle> mean_solar_time;
      ground_track.mean_solar_times_of_ascending_nodes_.emplace();
      for (auto ascending_node = ascending_nodes.Begin();
           ascending_node != ascending_nodes.End();
           ++ascending_node) {
        if (mean_solar_time.has_value()) {
          mean_solar_time = UnwindFrom(
              *mean_solar_time, MeanSolarTime(ascending_node, *mean_sun));
        } else {
          mean_solar_time =
              Mod(MeanSolarTime(ascending_node, *mean_sun), 2 * π * Radian);
        }
        ground_track.mean_solar_times_of_ascending_nodes_->Include(
            *mean_solar_time);
      }
    }
    if (nominal_recurrence.has_value()) {
      int n = 0;
      std::optional<Angle> initial_offset;
      std::optional<Angle> reduced_longitude;
      ground_track.reduced_longitudes_of_equator_crossings_of_ascending_passes_.emplace();
      for (auto ascending_node = ascending_nodes.Begin();
           ascending_node != ascending_nodes.End();
           ++ascending_node) {
        if (initial_offset.has_value()) {
          reduced_longitude = UnwindFrom(
              *reduced_longitude,
              PlanetocentricLongitude(ascending_node, primary) -
                  *initial_offset - n * nominal_recurrence->equatorial_shift());
        } else {
          auto const planetocentric_longitude =
              PlanetocentricLongitude(ascending_node, primary);
          reduced_longitude = Mod(planetocentric_longitude,
                                  nominal_recurrence->grid_interval());
          initial_offset = planetocentric_longitude - *reduced_longitude;
        }
        ground_track.reduced_longitudes_of_equator_crossings_of_ascending_passes_->Include(
            *reduced_longitude);
        ++n;
      }
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
