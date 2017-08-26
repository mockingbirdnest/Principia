
#include "ksp_plugin/interface.hpp"

#include "journal/method.hpp"
#include "journal/profiles.hpp"

namespace principia {
namespace interface {

std::future<void> const* principia__FutureCatchUpVessel(
    Plugin* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FutureCatchUpVessel> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(vessel_guid);
  auto future = plugin->CatchUpVessel(vessel_guid);
  return m.Return(future.release());
}

void principia__FutureWait(std::future<void> const** const future) {
  journal::Method<journal::FutureWait> m({future}, {future});
  TakeOwnership(future)->wait();
  return m.Return();
}

}  // namespace interface
}  // namespace principia
