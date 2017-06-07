
#include "ksp_plugin/planetarium.hpp"

#include "geometry/point.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::Position;
using quantities::Time;

Planetarium::Planetarium(
    std::vector<Sphere<Length, Barycentric>> const& spheres,
    Perspective<Barycentric, Camera, Length, OrthogonalMap> const& perspective,
    not_null<NavigationFrame*> const plotting_frame)
    : spheres_(spheres),
      perspective_(perspective),
      plotting_frame_(plotting_frame) {}

std::vector<RP2Point<Length, Camera>> Planetarium::PlotMethod0(
    DiscreteTrajectory<Barycentric> const& trajectory,
    Instant const& now) const {
  RigidMotion<Navigation, Barycentric> const inverse_rigid_motion_at_now =
      plotting_frame_->ToThisFrameAtTime(now).Inverse();
  std::vector<RP2Point<Length, Camera>> rp2_points;

  for (auto it = trajectory.Begin(); it != trajectory.End(); ++it) {
    AppendRP2PointIfNeeded(inverse_rigid_motion_at_now,
                           it.time(),
                           it.degrees_of_freedom(),
                           rp2_points);
  }

  return rp2_points;
}

std::vector<RP2Point<Length, Camera>> Planetarium::PlotMethod1(
    Trajectory<Barycentric> const& trajectory,
    Instant const& now,
    Length const& tolerance) const {
  RigidMotion<Navigation, Barycentric> const inverse_rigid_motion_at_now =
      plotting_frame_->ToThisFrameAtTime(now).Inverse();
  std::vector<RP2Point<Length, Camera>> rp2_points;

  Instant t = trajectory.t_min();
  Time Δt;
  while (t <= trajectory.t_max()) {
    DegreesOfFreedom<Barycentric> const barycentric_degrees_of_freedom =
        trajectory.EvaluateDegreesOfFreedom(t);
    AppendRP2PointIfNeeded(inverse_rigid_motion_at_now,
                           t,
                           barycentric_degrees_of_freedom,
                           rp2_points);

    // Don't pay any attention to the perspective, just adjust |Δt| to stay
    // within the |tolerance|.
    Δt = tolerance / barycentric_degrees_of_freedom.velocity().Norm();
    t += Δt;
  }

  return rp2_points;
}

void Planetarium::AppendRP2PointIfNeeded(
    RigidMotion<Navigation, Barycentric> const& inverse_rigid_motion_at_now,
    Instant const& t,
    DegreesOfFreedom<Barycentric> const& barycentric_degrees_of_freedom,
    std::vector<RP2Point<Length, Camera>>& rp2_points) const {
  RigidMotion<Barycentric, Navigation> const rigid_motion_at_t =
      plotting_frame_->ToThisFrameAtTime(t);
  RigidMotion<Barycentric, Barycentric> const rigid_motion_back_and_forth =
      inverse_rigid_motion_at_now * rigid_motion_at_t;
  DegreesOfFreedom<Barycentric> const plottable_degrees_of_freedom =
      rigid_motion_back_and_forth(barycentric_degrees_of_freedom);

  // TODO(phl): This is missing a precise determination of the time when the
  // trajectory become hidden.
  bool hidden = false;
  for (auto const& sphere : spheres_) {
    if (perspective_.IsHiddenBySphere(plottable_degrees_of_freedom.position(),
                                      sphere)) {
      hidden = true;
      break;
    }
  }
  if (!hidden) {
    RP2Point<Length, Camera> const rp2_point =
        perspective_(plottable_degrees_of_freedom.position());
    rp2_points.push_back(rp2_point);
  }
}

}  // namespace internal_planetarium
}  // namespace ksp_plugin
}  // namespace principia
