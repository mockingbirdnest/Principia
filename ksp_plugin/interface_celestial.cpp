#include "ksp_plugin/interface.hpp"

#include <string>
#include <vector>

#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "ksp_plugin/iterators.hpp"

namespace principia {
namespace interface {

using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_iterators;

Iterator* __cdecl principia__CelestialGetAllNames(
    Plugin const* const plugin) {
  journal::Method<journal::CelestialGetAllNames> m({plugin});
  CHECK(plugin != nullptr);
  auto const celestials = plugin->GetAllCelestials();
  std::vector<std::string> names;
  names.reserve(celestials.size());
  for (auto const& celestial : celestials) {
    names.push_back(celestial->body()->name());
  }
  return m.Return(new TypedIterator<std::vector<std::string>>(names, plugin));
}

// Calls `plugin->CelestialFromParent` with the arguments given.
// `plugin` must not be null.  No transfer of ownership.
QP __cdecl principia__CelestialFromParent(Plugin const* const plugin,
                                          int const celestial_index) {
  journal::Method<journal::CelestialFromParent> m({plugin, celestial_index});
  CHECK(plugin != nullptr);
  return m.Return(ToQP(plugin->CelestialFromParent(celestial_index)));
}

double __cdecl principia__CelestialInitialRotationInDegrees(
    Plugin const* const plugin,
    int const celestial_index) {
  journal::Method<journal::CelestialInitialRotationInDegrees> m(
      {plugin, celestial_index});
  CHECK(plugin != nullptr);
  return m.Return(plugin->CelestialInitialRotation(celestial_index) / Degree);
}

WXYZ __cdecl principia__CelestialRotation(Plugin const* const plugin,
                                          int const index) {
  journal::Method<journal::CelestialRotation> m({plugin, index});
  CHECK(plugin != nullptr);
  return m.Return(ToWXYZ(plugin->CelestialRotation(index).quaternion()));
}

double __cdecl principia__CelestialRotationPeriod(
    Plugin const* const plugin,
    int const celestial_index) {
  journal::Method<journal::CelestialRotationPeriod> m(
      {plugin, celestial_index});
  CHECK(plugin != nullptr);
  return m.Return(plugin->CelestialRotationPeriod(celestial_index) / Second);
}

WXYZ __cdecl principia__CelestialSphereRotation(Plugin const* const plugin) {
  journal::Method<journal::CelestialSphereRotation> m({plugin});
  CHECK(plugin != nullptr);
  return m.Return(ToWXYZ(plugin->CelestialSphereRotation().quaternion()));
}

QP __cdecl principia__CelestialWorldDegreesOfFreedom(Plugin const* const plugin,
                                                     int const index,
                                                     Origin const origin,
                                                     double const time) {
  journal::Method<journal::CelestialWorldDegreesOfFreedom> m(
      {plugin, index, origin, time});
  CHECK(plugin != nullptr);
  return m.Return(ToQP(
      plugin->CelestialWorldDegreesOfFreedom(
          index,
          plugin->BarycentricToWorld(
              origin.reference_part_is_unmoving,
              origin.reference_part_id,
              origin.reference_part_is_at_origin
                  ? std::nullopt
                  : std::make_optional(
                        FromXYZ<Position<World>>(
                            origin.main_body_centre_in_world))),
          FromGameTime(*plugin, time))));
}

}  // namespace interface
}  // namespace principia
