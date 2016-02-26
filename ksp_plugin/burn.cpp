
#include "ksp_plugin/burn.hpp"

namespace principia {
namespace ksp_plugin {

NavigationManœuvre MakeNavigationManœuvre(Burn burn, Mass const& initial_mass) {
  NavigationManœuvre manœuvre(burn.thrust,
                              initial_mass,
                              burn.specific_impulse,
                              NormalizeOrZero(burn.Δv),
                              std::move(burn.frame));
  manœuvre.set_initial_time(burn.initial_time);
  manœuvre.set_Δv(burn.Δv.Norm());
  return std::move(manœuvre);
}

}  // namespace ksp_plugin
}  // namespace principia
