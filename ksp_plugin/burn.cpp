
#include "ksp_plugin/burn.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_burn {

NavigationManœuvre MakeNavigationManœuvre(Burn burn, Mass const& initial_mass) {
  NavigationManœuvre::Intensity intensity;
  intensity.Δv = burn.Δv;
  NavigationManœuvre::Timing timing;
  timing.initial_time = burn.initial_time;
  return NavigationManœuvre(burn.thrust,
                            initial_mass,
                            burn.specific_impulse,
                            intensity,
                            timing,
                            std::move(burn.frame),
                            burn.is_inertially_fixed);
}

}  // namespace internal_burn
}  // namespace ksp_plugin
}  // namespace principia
