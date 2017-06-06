#include "planetarium.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

Planetarium::Planetarium(
    std::vector<Sphere<Length, Barycentric>> const& spheres,
    Perspective<Barycentric, Camera, Length, OrthogonalMap> const& perspective,
    not_null<NavigationFrame*> const plotting_frame)
    : spheres_(spheres),
      perspective_(perspective),
      plotting_frame_(plotting_frame) {}

std::vector<RP2Point<Length, Camera>> Planetarium::PlotMethod1(
    Trajectory<Barycentric> const& trajectory,
    Length const& tolerance) const {
  
}

}  // namespace internal_planetarium
}  // namespace ksp_plugin
}  // namespace principia
