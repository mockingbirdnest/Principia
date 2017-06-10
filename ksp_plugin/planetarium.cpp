
#include "ksp_plugin/planetarium.hpp"

#include "geometry/point.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::Position;
using quantities::Time;

Planetarium::Planetarium(
    std::vector<Sphere<Length, Barycentric>> const& spheres,
    Perspective<Navigation, Camera, Length, OrthogonalMap> const&
        perspective,
    not_null<NavigationFrame*> const plotting_frame)
    : spheres_(spheres),
      perspective_(perspective),
      plotting_frame_(plotting_frame) {}

std::vector<RP2Point<Length, Camera>> Planetarium::PlotMethod0(
    DiscreteTrajectory<Barycentric> const& trajectory,
    Instant const& now) const {
  auto const plottable_spheres = ComputePlottableSpheres(now);

  std::vector<RP2Point<Length, Camera>> rp2_points;
  for (auto it = trajectory.Begin(); it != trajectory.End(); ++it) {
    AppendRP2PointIfNeeded(it.time(),
                           it.degrees_of_freedom(),
                           plottable_spheres,
                           rp2_points);
  }
  return rp2_points;
}

std::vector<RP2Point<Length, Camera>> Planetarium::PlotMethod1(
    Trajectory<Barycentric> const& trajectory,
    Instant const& now,
    Length const& tolerance) const {
  auto const plottable_spheres = ComputePlottableSpheres(now);

  std::vector<RP2Point<Length, Camera>> rp2_points;
  Instant t = trajectory.t_min();
  Time Δt;
  while (t <= trajectory.t_max()) {
    DegreesOfFreedom<Barycentric> const barycentric_degrees_of_freedom =
        trajectory.EvaluateDegreesOfFreedom(t);
    AppendRP2PointIfNeeded(t,
                           barycentric_degrees_of_freedom,
                           plottable_spheres,
                           rp2_points);

    // Don't pay any attention to the perspective, just adjust |Δt| to stay
    // within the |tolerance|.
    Δt = tolerance / barycentric_degrees_of_freedom.velocity().Norm();
    t += Δt;
  }
  return rp2_points;
}

std::vector<Sphere<Length, Navigation>> Planetarium::ComputePlottableSpheres(
    Instant const& now) const {
  RigidMotion<Barycentric, Navigation> const rigid_motion_at_now =
      plotting_frame_->ToThisFrameAtTime(now);

  std::vector<Sphere<Length, Navigation>> plottable_spheres;
  for (auto const& barycentric_sphere : spheres_) {
    plottable_spheres.emplace_back(
        rigid_motion_at_now.rigid_transformation()(barycentric_sphere.centre()),
        barycentric_sphere.radius());
  }
  return plottable_spheres;
}

void Planetarium::AppendRP2PointIfNeeded(
    Instant const& t,
    DegreesOfFreedom<Barycentric> const& barycentric_degrees_of_freedom,
    std::vector<Sphere<Length, Navigation>> const& plottable_spheres,
    std::vector<RP2Point<Length, Camera>>& rp2_points) const {
  RigidMotion<Barycentric, Navigation> const rigid_motion_at_t =
      plotting_frame_->ToThisFrameAtTime(t);
  DegreesOfFreedom<Navigation> const plottable_degrees_of_freedom =
      rigid_motion_at_t(barycentric_degrees_of_freedom);

  // TODO(phl): This happily draws behind the camera.
  // TODO(phl): This is missing a precise determination of the time when the
  // trajectory become hidden.
  bool hidden = false;
  for (auto const& plottable_sphere : plottable_spheres) {
    if (perspective_.IsHiddenBySphere(plottable_degrees_of_freedom.position(),
                                      plottable_sphere)) {
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
