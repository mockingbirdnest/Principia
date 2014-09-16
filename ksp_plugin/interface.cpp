#include "ksp_plugin/interface.hpp"

#include <string>

using principia::si::Degree;
using principia::si::Metre;

namespace principia {
namespace ksp_plugin {

Plugin* CreatePlugin(Instant const initial_time, int const sun_index,
                     GravitationalParameter const sun_gravitational_parameter,
                     double const planetarium_rotation_in_degrees) {
  return new Plugin(initial_time, sun_index, sun_gravitational_parameter,
                    planetarium_rotation_in_degrees * Degree);
}

void DestroyPlugin(Plugin* plugin) {
  delete plugin;
  plugin = nullptr;
}

void InsertCelestial(Plugin* plugin, int const index,
                     GravitationalParameter const gravitational_parameter,
                     int const parent,
                     Displacement<AliceSun> const from_parent_position,
                     Velocity<AliceSun> const from_parent_velocity) {
  plugin->InsertCelestial(index, gravitational_parameter, parent,
                          from_parent_position, from_parent_velocity);
}

void UpdateCelestialHierarchy(Plugin* plugin, int const index,
                              int const parent) {
  plugin->UpdateCelestialHierarchy(index, parent);
}

void InsertOrKeepVessel(Plugin* plugin, char const* guid, int const parent) {
  plugin->InsertOrKeepVessel(guid, parent);
}

void SetVesselStateOffset(Plugin* plugin, char const* guid,
                          Displacement<AliceSun> const from_parent_position,
                          Velocity<AliceSun> const from_parent_velocity) {
  plugin->SetVesselStateOffset(guid, from_parent_position,
                               from_parent_velocity);
}

XYZ VesselDisplacementFromParent(Plugin* plugin, char const* guid) {
  R3Element<Length> const result =
      plugin->VesselDisplacementFromParent(guid).coordinates();
  return {result.x / Metre, result.y / Metre, result.z / Metre};
}

XYZ VesselParentRelativeVelocity(Plugin* plugin, char const* guid) {
  R3Element<Speed> const result =
      plugin->VesselParentRelativeVelocity(guid).coordinates();
  return {result.x / (Metre / Second),
          result.y / (Metre / Second),
          result.z / (Metre / Second)};
}

XYZ CelestialDisplacementFromParent(Plugin* plugin, int const index) {
  R3Element<Length> const result =
      plugin->CelestialDisplacementFromParent(index).coordinates();
  return {result.x / Metre, result.y / Metre, result.z / Metre};;
}

XYZ CelestialParentRelativeVelocity(Plugin* plugin, int const index) {
  R3Element<Speed> const result =
      plugin->CelestialParentRelativeVelocity(index).coordinates();
  return {result.x / (Metre / Second),
          result.y / (Metre / Second),
          result.z / (Metre / Second)};
}

char const* SayHello() {
  return "Hello from native C++!";
}

}  // namespace ksp_plugin
}  // namespace principia
