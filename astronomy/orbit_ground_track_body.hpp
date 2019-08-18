#pragma once

#include "astronomy/orbit_ground_track.hpp"

#include "physics/apsides.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_ground_track {

using geometry::Vector;
using physics::ComputeNodes;
using quantities::Mod;
using quantities::UnwindFrom;
using quantities::si::Radian;

template<typename PrimaryCentred, typename Inertial>
OrbitGroundTrack OrbitGroundTrack::ForTrajectory(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    RotatingBody<Inertial> const& primary,
    std::optional<OrbitRecurrence> const& nominal_recurrence) {
  DiscreteTrajectory<PrimaryCentred> ascending_nodes;
  DiscreteTrajectory<PrimaryCentred> descending_nodes;
  ComputeNodes(trajectory.Begin(),
               trajectory.End(),
               Vector<double, PrimaryCentred>({0, 0, 1}),
               /*max_points=*/std::numeric_limits<int>::max(),
               ascending_nodes,
               descending_nodes);
  OrbitGroundTrack ground_track;
  if (nominal_recurrence.has_value() && ascending_nodes.Size() > 0) {
    int n = 0;
    std::optional<Angle> initial_offset;
    std::optional<Angle> reduced_longitude;
    ground_track.reduced_longitude_of_equator_crossing_.emplace();
    for (auto ascending_node = ascending_nodes.Begin();
         ascending_node != ascending_nodes.End();
         ++ascending_node) {
      Angle const celestial_longitude =
          (ascending_node.degrees_of_freedom().position() -
           PrimaryCentred::origin)
              .coordinates()
              .ToSpherical()
              .longitude;
      Angle const planetocentric_longitude =
          celestial_longitude - primary.AngleAt(ascending_node.time()) -
          π / 2 * Radian;
      if (initial_offset.has_value()) {
        reduced_longitude =
            UnwindFrom(*reduced_longitude,
                       planetocentric_longitude - *initial_offset -
                           n * nominal_recurrence->equatorial_shift());
      } else {
        reduced_longitude =
            Mod(planetocentric_longitude, nominal_recurrence->grid_interval());
        initial_offset = planetocentric_longitude - *reduced_longitude;
      }
      ground_track.reduced_longitude_of_equator_crossing_->Include(
          *reduced_longitude);
      ++n;
    }
  }
  return ground_track;
}

inline std::optional<Interval<Angle>> const&
OrbitGroundTrack::reduced_longitude_of_equator_crossing() const {
  return reduced_longitude_of_equator_crossing_;
}

}  // namespace internal_orbit_ground_track
}  // namespace astronomy
}  // namespace principia
