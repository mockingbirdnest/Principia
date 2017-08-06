
#include "ksp_plugin_test/mock_plugin.hpp"

#include <vector>

#include "gmock/gmock.h"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

MockPlugin::MockPlugin() : Plugin("JD2451545", "JD2451545", Angle()) {}

not_null<std::unique_ptr<Planetarium>> MockPlugin::NewPlanetarium(
    Planetarium::Parameters const& parameters,
    Perspective<Navigation, Camera> const& perspective) const {
  std::unique_ptr<Planetarium> planetarium;
  FillPlanetarium(parameters, perspective, &planetarium);
  return std::move(planetarium);
}

not_null<std::unique_ptr<NavigationFrame>>
MockPlugin::NewBodyCentredNonRotatingNavigationFrame(
    Index const reference_body_index) const {
  std::unique_ptr<NavigationFrame> navigation_frame;
  FillBodyCentredNonRotatingNavigationFrame(reference_body_index,
                                            &navigation_frame);
  return std::move(navigation_frame);
}

not_null<std::unique_ptr<NavigationFrame>>
MockPlugin::NewBarycentricRotatingNavigationFrame(
    Index const primary_index,
    Index const secondary_index) const {
  std::unique_ptr<NavigationFrame> navigation_frame;
  FillBarycentricRotatingNavigationFrame(primary_index,
                                         secondary_index,
                                         &navigation_frame);
  return std::move(navigation_frame);
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
