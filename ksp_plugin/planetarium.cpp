
#include "ksp_plugin/planetarium.hpp"

#include <vector>

#include "geometry/point.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::Position;
using geometry::RP2Line;
using quantities::Pow;
using quantities::Sin;
using quantities::Tan;
using quantities::Time;

Planetarium::Parameters::Parameters(double const sphere_radius_multiplier,
                                    Angle const& angular_resolution,
                                    Angle const& field_of_view)
    : sphere_radius_multiplier_(sphere_radius_multiplier),
      sin²_angular_resolution_(Pow<2>(Sin(angular_resolution))),
      tan_field_of_view_(Tan(field_of_view)) {}

Planetarium::Planetarium(
    Parameters const& parameters,
    Perspective<Navigation, Camera, Length, OrthogonalMap> const& perspective,
    not_null<Ephemeris<Barycentric> const*> const ephemeris,
    not_null<NavigationFrame const*> const plotting_frame)
    : parameters_(parameters),
      perspective_(perspective),
      ephemeris_(ephemeris),
      plotting_frame_(plotting_frame) {}

RP2Lines<Length, Camera> Planetarium::PlotMethod0(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now) const {
  auto const plottable_spheres = ComputePlottableSpheres(now);
  auto const plottable_segments = ComputePlottableSegments(plottable_spheres,
                                                           begin, end);

  auto const field_of_view_radius² =
      perspective_.focal() * perspective_.focal() *
      parameters_.tan_field_of_view_ * parameters_.tan_field_of_view_;
  RP2Lines<Length, Camera> rp2_lines;
  for (auto const& plottable_segment : plottable_segments) {
    // Apply the projection to the current plottable segment.
    auto const rp2_first = perspective_(plottable_segment.first);
    auto const rp2_second = perspective_(plottable_segment.second);

    // If the segment is entirely outside the field of view, ignore it.
    Length const x1 = rp2_first.x();
    Length const y1 = rp2_first.y();
    Length const x2 = rp2_second.x();
    Length const y2 = rp2_second.y();
    if (x1 * x1 + y1 * y1 > field_of_view_radius² &&
        x2 * x2 + y2 * y2 > field_of_view_radius²) {
      continue;
    }

    // Create a new ℝP² line when two segments are not consecutive.
    if (rp2_lines.empty() || rp2_lines.back().back() != rp2_first) {
      RP2Line<Length, Camera> const rp2_line = {rp2_first, rp2_second};
      rp2_lines.push_back(rp2_line);
    } else {
      rp2_lines.back().push_back(rp2_second);
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
    Sphere<Length, Navigation> plottable_sphere(
        rigid_motion_at_now.rigid_transformation()(centre_in_barycentric),
        parameters_.sphere_radius_multiplier_ * mean_radius);
    // If the sphere is seen under an angle that is very small it doesn't
    // participate in hiding.
    if (perspective_.SphereSin²HalfAngle(plottable_sphere) <
        parameters_.sin²_angular_resolution_) {
      plottable_spheres.emplace_back(std::move(plottable_sphere));
    }
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

    // Find the part of the segment that is behind the focal plane.  We don't
    // care about things that are in front of the focal plane.
    const Segment<Displacement<Navigation>> segment = {p1, p2};
    auto const segment_behind_focal_plane =
        perspective_.SegmentBehindFocalPlane(segment);
    if (segment_behind_focal_plane) {
      // Find the part(s) of the segment that are not hidden by spheres.  These
      // are the ones we want to plot.
      auto segments = perspective_.VisibleSegments(*segment_behind_focal_plane,
                                                   plottable_spheres);
      std::move(segments.begin(),
                segments.end(),
                std::back_inserter(all_segments));
    }

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
