
#include "ksp_plugin/burn.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_burn {

NavigationManœuvre MakeNavigationManœuvre(Burn burn, Mass const& initial_mass) {
  NavigationManœuvre manœuvre(burn.thrust,
                              initial_mass,
                              burn.specific_impulse,
                              NormalizeOrZero(burn.Δv),
                              std::move(burn.frame),
                              burn.is_inertially_fixed);
  manœuvre.set_initial_time(burn.initial_time);
  manœuvre.set_Δv(burn.Δv.Norm());
  return manœuvre;
}

}  // namespace internal_burn
}  // namespace ksp_plugin
}  // namespace principia
