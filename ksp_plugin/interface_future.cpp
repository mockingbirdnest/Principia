
#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <chrono>

#include "journal/method.hpp"
#include "journal/profiles.hpp"

namespace principia {
namespace interface {

PileUpFuture const* principia__FutureCatchUpVessel(
    Plugin* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FutureCatchUpVessel> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(vessel_guid);
  auto future = plugin->CatchUpVessel(vessel_guid);
  return m.Return(future.release());
}

void principia__FutureWait(PileUpFuture const** const future) {
  journal::Method<journal::FutureWait> m({future}, {future});
  auto const owned_future = TakeOwnership(future);
  owned_future->future.wait();
  return m.Return();
}

}  // namespace interface
}  // namespace principia
