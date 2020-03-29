
#include "ksp_plugin/planetarium.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

#include "geometry/point.hpp"
#include "physics/massive_body.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::Position;
using geometry::RP2Line;
using geometry::Sign;
using geometry::Velocity;
using physics::MassiveBody;
using quantities::Pow;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Tan;
using quantities::Time;

namespace {
constexpr int max_plot_method_2_steps = 10'000;
}  // namespace

Planetarium::Parameters::Parameters(double const sphere_radius_multiplier,
                                    Angle const& angular_resolution,
                                    Angle const& field_of_view)
    : sphere_radius_multiplier_(sphere_radius_multiplier),
      sin²_angular_resolution_(Pow<2>(Sin(angular_resolution))),
      tan_angular_resolution_(Tan(angular_resolution)),
      tan_field_of_view_(Tan(field_of_view)) {}

Planetarium::Planetarium(
    Parameters const& parameters,
    Perspective<Navigation, Camera> perspective,
    not_null<Ephemeris<Barycentric> const*> const ephemeris,
    not_null<NavigationFrame const*> const plotting_frame)
    : parameters_(parameters),
      perspective_(std::move(perspective)),
      ephemeris_(ephemeris),
      plotting_frame_(plotting_frame) {}

RP2Lines<Length, Camera> Planetarium::PlotMethod0(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now,
    bool const /*reverse*/) const {
  auto const plottable_begin =
      begin.trajectory()->LowerBound(plotting_frame_->t_min());
  auto const plottable_end =
      begin.trajectory()->LowerBound(plotting_frame_->t_max());
  auto const plottable_spheres = ComputePlottableSpheres(now);
  auto const plottable_segments = ComputePlottableSegments(plottable_spheres,
                                                           plottable_begin,
                                                           plottable_end);

  auto const field_of_view_radius² =
      perspective_.focal() * perspective_.focal() *
      parameters_.tan_field_of_view_ * parameters_.tan_field_of_view_;
  std::optional<Position<Navigation>> previous_position;
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

    // Create a new ℝP² line when two segments are not consecutive.  Don't
    // compare ℝP² points for equality, that's expensive.
    bool const are_consecutive =
        previous_position == plottable_segment.first;
    previous_position = plottable_segment.second;

    if (are_consecutive) {
      rp2_lines.back().push_back(rp2_second);
    } else {
      RP2Line<Length, Camera> const rp2_line = {rp2_first, rp2_second};
      rp2_lines.push_back(rp2_line);
    }
  }
  return rp2_lines;
}

RP2Lines<Length, Camera> Planetarium::PlotMethod1(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now,
    bool const reverse) const {
  Length const focal_plane_tolerance =
      perspective_.focal() * parameters_.tan_angular_resolution_;
  auto const focal_plane_tolerance² =
      focal_plane_tolerance * focal_plane_tolerance;

  auto const rp2_lines = PlotMethod0(begin, end, now, reverse);

  RP2Lines<Length, Camera> new_rp2_lines;
  for (auto const& rp2_line : rp2_lines) {
    RP2Line<Length, Camera> new_rp2_line;
    std::optional<RP2Point<Length, Camera>> start_rp2_point;
    for (int i = 0; i < rp2_line.size(); ++i) {
      RP2Point<Length, Camera> const& rp2_point = rp2_line[i];
      if (i == 0) {
        new_rp2_line.push_back(rp2_point);
        start_rp2_point = rp2_point;
      } else if (Pow<2>(rp2_point.x() - start_rp2_point->x()) +
                 Pow<2>(rp2_point.y() - start_rp2_point->y()) >
                     focal_plane_tolerance²) {
        // TODO(phl): This creates a segment if the tolerance is exceeded.  It
        // should probably create a segment that stays just below the tolerance.
        new_rp2_line.push_back(rp2_point);
        start_rp2_point = rp2_point;
      } else if (i == rp2_line.size() - 1) {
        new_rp2_line.push_back(rp2_point);
      }
    }
    new_rp2_lines.push_back(std::move(new_rp2_line));
  }
  return new_rp2_lines;
}

RP2Lines<Length, Camera> Planetarium::PlotMethod2(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now,
    bool const reverse) const {
  if (begin == end) {
    return {};
  }
  auto last = end;
  --last;
  auto const& trajectory = *begin.trajectory();
  auto const begin_time = std::max(begin->time, plotting_frame_->t_min());
  auto const last_time = std::min(last->time, plotting_frame_->t_max());
  return PlotMethod2(trajectory, begin_time, last_time, now, reverse);
}

RP2Lines<Length, Camera> Planetarium::PlotMethod2(
    Trajectory<Barycentric> const& trajectory,
    Instant const& first_time,
    Instant const& last_time,
    Instant const& now,
    bool const reverse) const {
  RP2Lines<Length, Camera> lines;
  auto const plottable_spheres = ComputePlottableSpheres(now);
  double const tan²_angular_resolution =
      Pow<2>(parameters_.tan_angular_resolution_);
  auto const final_time = reverse ? first_time : last_time;
  auto previous_time = reverse ? last_time : first_time;

  Sign const direction = reverse ? Sign::Negative() : Sign::Positive();
  if (direction * (final_time - previous_time) <= Time{}) {
    return lines;
  }
  RigidMotion<Barycentric, Navigation> to_plotting_frame_at_t =
      plotting_frame_->ToThisFrameAtTime(previous_time);
  DegreesOfFreedom<Navigation> const initial_degrees_of_freedom =
      to_plotting_frame_at_t(
          trajectory.EvaluateDegreesOfFreedom(previous_time));
  Position<Navigation> previous_position =
      initial_degrees_of_freedom.position();
  Velocity<Navigation> previous_velocity =
      initial_degrees_of_freedom.velocity();
  Time Δt = final_time - previous_time;

  Instant t;
  double estimated_tan²_error;
  std::optional<DegreesOfFreedom<Barycentric>>
      degrees_of_freedom_in_barycentric;
  Position<Navigation> position;

  std::optional<Position<Navigation>> last_endpoint;

  int steps_accepted = 0;

  goto estimate_tan²_error;

  while (steps_accepted < max_plot_method_2_steps &&
         direction * (previous_time - final_time) < Time{}) {
    do {
      // One square root because we have squared errors, another one because the
      // errors are quadratic in time (in other words, two square roots because
      // the squared errors are quartic in time).
      // A safety factor prevents catastrophic retries.
      Δt *= 0.9 * Sqrt(Sqrt(tan²_angular_resolution / estimated_tan²_error));
    estimate_tan²_error:
      t = previous_time + Δt;
      if (direction * (t - final_time) > Time{}) {
        t = final_time;
        Δt = t - previous_time;
      }
      Position<Navigation> const extrapolated_position =
          previous_position + previous_velocity * Δt;
      to_plotting_frame_at_t = plotting_frame_->ToThisFrameAtTime(t);
      degrees_of_freedom_in_barycentric =
          trajectory.EvaluateDegreesOfFreedom(t);
      position = to_plotting_frame_at_t.rigid_transformation()(
                     degrees_of_freedom_in_barycentric->position());

      // The quadratic term of the error between the linear interpolation and
      // the actual function is maximized halfway through the segment, so it is
      // 1/2 (Δt/2)² f″(t-Δt) = (1/2 Δt² f″(t-Δt)) / 4; the squared error is
      // thus (1/2 Δt² f″(t-Δt))² / 16.
      estimated_tan²_error =
          perspective_.Tan²AngularDistance(extrapolated_position, position) /
          16;
    } while (estimated_tan²_error > tan²_angular_resolution);
    ++steps_accepted;

    // TODO(egg): also limit to field of view.
    auto const segment_behind_focal_plane =
        perspective_.SegmentBehindFocalPlane(
            Segment<Navigation>(previous_position, position));

    previous_time = t;
    previous_position = position;
    previous_velocity =
        to_plotting_frame_at_t(*degrees_of_freedom_in_barycentric).velocity();

    if (!segment_behind_focal_plane) {
      continue;
    }

    auto const visible_segments = perspective_.VisibleSegments(
                                      *segment_behind_focal_plane,
                                      plottable_spheres);
    for (auto const& segment : visible_segments) {
      if (last_endpoint != segment.first) {
        lines.emplace_back();
        lines.back().push_back(perspective_(segment.first));
      }
      lines.back().push_back(perspective_(segment.second));
      last_endpoint = segment.second;
    }
  }
  return lines;
}

std::vector<Sphere<Navigation>> Planetarium::ComputePlottableSpheres(
    Instant const& now) const {
  RigidMotion<Barycentric, Navigation> const rigid_motion_at_now =
      plotting_frame_->ToThisFrameAtTime(now);
  std::vector<Sphere<Navigation>> plottable_spheres;

  auto const& bodies = ephemeris_->bodies();
  for (not_null<MassiveBody const*> const body : bodies) {
    auto const trajectory = ephemeris_->trajectory(body);
    Length const mean_radius = body->mean_radius();
    Position<Barycentric> const centre_in_barycentric =
        trajectory->EvaluatePosition(now);
    Sphere<Navigation> plottable_sphere(
        rigid_motion_at_now.rigid_transformation()(centre_in_barycentric),
        parameters_.sphere_radius_multiplier_ * mean_radius);
    // If the sphere is seen under an angle that is very small it doesn't
    // participate in hiding.
    if (perspective_.SphereSin²HalfAngle(plottable_sphere) >
        parameters_.sin²_angular_resolution_) {
      plottable_spheres.emplace_back(std::move(plottable_sphere));
    }
  }
  return plottable_spheres;
}

Segments<Navigation> Planetarium::ComputePlottableSegments(
    const std::vector<Sphere<Navigation>>& plottable_spheres,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end) const {
  Segments<Navigation> all_segments;
  if (begin == end) {
    return all_segments;
  }
  auto it1 = begin;
  Instant t1 = it1->time;
  RigidMotion<Barycentric, Navigation> rigid_motion_at_t1 =
      plotting_frame_->ToThisFrameAtTime(t1);
  Position<Navigation> p1 =
      rigid_motion_at_t1(it1->degrees_of_freedom).position();

  auto it2 = it1;
  while (++it2 != end) {
    // Processing one segment of the trajectory.
    Instant const t2 = it2->time;

    // Transform the degrees of freedom to the plotting frame.
    RigidMotion<Barycentric, Navigation> const rigid_motion_at_t2 =
        plotting_frame_->ToThisFrameAtTime(t2);
    Position<Navigation> const p2 =
        rigid_motion_at_t2(it2->degrees_of_freedom).position();

    // Find the part of the segment that is behind the focal plane.  We don't
    // care about things that are in front of the focal plane.
    const Segment<Navigation> segment = {p1, p2};
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
