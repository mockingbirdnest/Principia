
#include "ksp_plugin/planetarium.hpp"

#include <vector>

#include "geometry/point.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::Position;
using quantities::Time;

Planetarium::Parameters::Parameters(double const sphere_radius_multiplier)
    : sphere_radius_multiplier_(sphere_radius_multiplier) {}

Planetarium::Planetarium(
    Parameters const& parameters,
    Perspective<Navigation, Camera, Length, OrthogonalMap> const& perspective,
    not_null<Ephemeris<Barycentric> const*> const ephemeris,
    not_null<NavigationFrame*> const plotting_frame)
    : parameters_(parameters),
      perspective_(perspective),
      ephemeris_(ephemeris),
      plotting_frame_(plotting_frame) {}

std::vector<RP2Line<Length, Camera>> Planetarium::PlotMethod0(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now) const {
  auto const plottable_spheres = ComputePlottableSpheres(now);
  auto const plottable_segments = ComputePlottableSegments(plottable_spheres,
                                                           begin, end);

  std::vector<RP2Line<Length, Camera>> rp2_lines;
  for (auto const& plottable_segment : plottable_segments) {
    // Apply the projection to the current plottable segment.
    auto const rp2_first = perspective_(plottable_segment.first);
    auto const rp2_second = perspective_(plottable_segment.second);

    // Ignore any segment that goes behind the camera.
    if (!rp2_first || !rp2_second) {
      continue;
    }

    // Create a new ℝP² line when two segments are not consecutive.
    if (rp2_lines.empty() || rp2_lines.back().back() != *rp2_first) {
      RP2Line<Length, Camera> const rp2_line = {*rp2_first, *rp2_second};
      rp2_lines.push_back(rp2_line);
    } else {
      rp2_lines.back().push_back(*rp2_second);
    }
  }
  return rp2_lines;
}

std::vector<Sphere<Length, Navigation>> Planetarium::ComputePlottableSpheres(
    Instant const& now) const {
  RigidMotion<Barycentric, Navigation> const rigid_motion_at_now =
      plotting_frame_->ToThisFrameAtTime(now);
  std::vector<Sphere<Length, Navigation>> plottable_spheres;

  auto const& bodies = ephemeris_->bodies();
  for (auto const body : bodies) {
    auto const trajectory = ephemeris_->trajectory(body);
    Length const mean_radius = body->mean_radius();
    Position<Barycentric> const centre_in_barycentric =
        trajectory->EvaluatePosition(now);
    // TODO(phl): Don't create a plottable sphere if the body is very far from
    // the camera.  What should the criteria be?
    plottable_spheres.emplace_back(
        rigid_motion_at_now.rigid_transformation()(centre_in_barycentric),
        parameters_.sphere_radius_multiplier_ * mean_radius);
  }
  return plottable_spheres;
}

std::vector<Segment<Displacement<Navigation>>>
Planetarium::ComputePlottableSegments(
    const std::vector<Sphere<Length, Navigation>>& plottable_spheres,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end) const {
  std::vector<Segment<Displacement<Navigation>>> all_segments;
  if (begin == end) {
    return all_segments;
  }
  auto it1 = begin;
  Instant t1 = it1.time();
  RigidMotion<Barycentric, Navigation> rigid_motion_at_t1 =
      plotting_frame_->ToThisFrameAtTime(t1);
  Position<Navigation> p1 =
      rigid_motion_at_t1(it1.degrees_of_freedom()).position();

  auto it2 = it1;
  while (++it2 != end) {
    // Processing one segment of the trajectory.
    Instant const t2 = it2.time();

    // Transform the degrees of freedom to the plotting frame.
    RigidMotion<Barycentric, Navigation> const rigid_motion_at_t2 =
        plotting_frame_->ToThisFrameAtTime(t2);
    Position<Navigation> const p2 =
        rigid_motion_at_t2(it2.degrees_of_freedom()).position();

    // Use the perspective to compute the visible segments in the Navigation
    // frame.
    const Segment<Displacement<Navigation>> segment = {p1, p2};
    auto segments = perspective_.VisibleSegments(segment, plottable_spheres);
    std::move(segments.begin(),
              segments.end(),
              std::back_inserter(all_segments));

    it1 = it2;
    t1 = t2;
    rigid_motion_at_t1 = rigid_motion_at_t2;
    p1 = p2;
  }

  return all_segments;
}

}  // namespace internal_planetarium
}  // namespace ksp_plugin
}  // namespace principia
