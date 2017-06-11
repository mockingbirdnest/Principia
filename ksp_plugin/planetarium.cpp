
#include "ksp_plugin/planetarium.hpp"

#include <vector>

#include "geometry/point.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::Position;
using quantities::Time;

namespace {

constexpr double sphere_radius_multiplier = 1.05;

}  // namespace

Planetarium::Planetarium(
    Perspective<Navigation, Camera, Length, OrthogonalMap> const& perspective,
    not_null<Ephemeris<Barycentric> const*> const ephemeris,
    not_null<NavigationFrame*> const plotting_frame)
    : perspective_(perspective),
      ephemeris_(ephemeris),
      plotting_frame_(plotting_frame) {}

std::vector<RP2Point<Length, Camera>> Planetarium::PlotMethod0(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now) const {
  auto const plottable_spheres = ComputePlottableSpheres(now);

  std::vector<RP2Point<Length, Camera>> rp2_points;
  for (auto it = begin; it != end; ++it) {
    AppendRP2PointIfNeeded(it.time(),
                           it.degrees_of_freedom(),
                           plottable_spheres,
                           rp2_points);
  }
  return rp2_points;
}

std::vector<RP2Point<Length, Camera>> Planetarium::PlotMethod1(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now,
    Length const& tolerance) const {
  auto const plottable_spheres = ComputePlottableSpheres(now);

  auto const& trajectory = *begin.trajectory();
  std::vector<RP2Point<Length, Camera>> rp2_points;
  Instant t = begin.time();
  Time Δt;
  while (t < end.time()) {
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

  auto const& bodies = ephemeris_->bodies();
  for (auto const body : bodies) {
    auto const trajectory = ephemeris_->trajectory(body);
    Length const mean_radius = body->mean_radius();
    Position<Barycentric> const centre_in_barycentric =
        trajectory->EvaluatePosition(now);
    plottable_spheres.emplace_back(
        rigid_motion_at_now.rigid_transformation()(centre_in_barycentric),
        sphere_radius_multiplier * mean_radius);
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
