#pragma once

#include "ksp_plugin/plugin.hpp"

// DLL-exported functions for interfacing with Platform Invocation Services.

#if defined(_WIN32)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __attribute__((visibility("default")))
#endif

namespace principia {
namespace ksp_plugin {

struct XYZ { double x, y, z; };

// Returns a pointer to a plugin constructed with the arguments given.
// |planetarium_rotation_in_degrees| is converted to an |Angle| (whose
// underlying value is in radians), the other parameters are passed unmodified.
extern "C" DLLEXPORT
Plugin* CreatePlugin(Instant const initial_time, int const sun_index,
                     GravitationalParameter const sun_gravitational_parameter,
                     double const planetarium_rotation_in_degrees);

// Deletes |plugin|.
extern "C" DLLEXPORT
void DestroyPlugin(Plugin* plugin);

// Calls |plugin->InsertCelestial| with the arguments given. The arguments are
// passed unmodified.
extern "C" DLLEXPORT
void InsertCelestial(Plugin* plugin, int const index,
                     GravitationalParameter const gravitational_parameter,
                     int const parent,
                     Displacement<AliceSun> const from_parent_position,
                     Velocity<AliceSun> const from_parent_velocity);

// Calls |plugin->InsertCelestial| with the arguments given. The arguments are
// passed unmodified.
extern "C" DLLEXPORT
void UpdateCelestialHierarchy(Plugin* plugin, int const index,
                              int const parent);

// Calls |plugin->InsertOrKeepVessel| with the arguments given. |guid| is
// converted to |std::string|, |parent| is passed unmodified.
extern "C" DLLEXPORT
void InsertOrKeepVessel(Plugin* plugin, char const* guid, int const parent);

// Calls |plugin->SetVesselStateOffset| with the arguments given. |guid| is
// converted to |std::string|, the other parameters are passed unmodified.
extern "C" DLLEXPORT
void SetVesselStateOffset(Plugin* plugin, char const* guid,
                          Displacement<AliceSun> const from_parent_position,
                          Velocity<AliceSun> const from_parent_velocity);

// Calls |plugin->VesselDisplacementFromParent| with the arguments given. |guid|
// is converted to |std::string|, the other parameters are passed unmodified.
extern "C" DLLEXPORT
XYZ VesselDisplacementFromParent(Plugin* plugin,
                                 char const* guid);

// Calls |plugin->VesselParentRelativeVelocity| with the arguments given. |guid|
// is converted to |std::string|, the other parameters are passed unmodified.
extern "C" DLLEXPORT
XYZ VesselParentRelativeVelocity(Plugin* plugin,
                                 char const* guid);

// Calls |plugin->CelestialDisplacementFromParent| with the arguments given. The
// parameters are passed unmodified.
extern "C" DLLEXPORT
XYZ CelestialDisplacementFromParent(Plugin* plugin, int const index);

// Calls |plugin->CelestialParentRelativeVelocity| with the arguments given. The
// parameters are passed unmodified.
extern "C" DLLEXPORT
XYZ CelestialParentRelativeVelocity(Plugin* plugin, int const index);

extern "C" DLLEXPORT
char const* SayHello();

}  // namespace ksp_plugin
}  // namespace principia
