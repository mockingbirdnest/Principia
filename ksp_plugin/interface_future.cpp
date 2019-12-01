
#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <chrono>

#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/iterators.hpp"

namespace principia {
namespace interface {

using ksp_plugin::TypedIterator;
using ksp_plugin::VesselSet;

PileUpFuture* __cdecl principia__FutureCatchUpVessel(
    Plugin* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FutureCatchUpVessel> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(vessel_guid);
  auto future = plugin->CatchUpVessel(vessel_guid);
  return m.Return(future.release());
}

void __cdecl principia__FutureWaitForVesselToCatchUp(
    Plugin* const plugin,
    PileUpFuture** const future,
    Iterator** const collided_vessels) {
  journal::Method<journal::FutureWaitForVesselToCatchUp> m(
      {plugin, future}, {future, collided_vessels});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(future);
  auto const owned_future = TakeOwnership(future);
  VesselSet collided_vessel_set;
  plugin->WaitForVesselToCatchUp(*owned_future, collided_vessel_set);
  *collided_vessels =
      new TypedIterator<VesselSet>(std::move(collided_vessel_set));
  return m.Return();
}

}  // namespace interface
}  // namespace principia
