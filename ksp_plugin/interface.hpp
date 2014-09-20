#pragma once

#include <type_traits>

#include "ksp_plugin/plugin.hpp"

// DLL-exported functions for interfacing with Platform Invocation Services.

#if defined(_WIN32)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __attribute__((visibility("default")))
#endif

namespace principia {
namespace ksp_plugin {

extern "C"
struct XYZ { double x, y, z; };
static_assert(std::is_standard_layout<XYZ>::value,
              "XYZ is used for interfacing");

// Sets stderr to log INFO, and redirects stderr, which Unity does not log, to
// "<KSP directory>/stderr.log". This provides an easily accessible file
// containing a sufficiently verbose log of the latest session, instead of
// requiring users to dig in the archive of all past logs at all severities.
// This archive is written to
// "<KSP directory>/glog/Principia/<SEVERITY>.<date>-<time>.<pid>", 
// where date and time are in ISO 8601 basic format.
// TODO(egg): libglog should really be statically linked, what happens if two
// plugins use glog?
extern "C" DLLEXPORT
void InitGoogleLogging();

// Exports |LOG(SEVERITY) << message| for fast logging from the C# adapter.
extern "C" DLLEXPORT
void LOGINFO(char const* message);
extern "C" DLLEXPORT
void LOGWARNING(char const* message);
extern "C" DLLEXPORT
void LOGERROR(char const* message);
extern "C" DLLEXPORT
void LOGFATAL(char const* message);

// Returns a pointer to a plugin constructed with the arguments given.
extern "C" DLLEXPORT
Plugin* CreatePlugin(double const initial_time, int const sun_index,
                     double const sun_gravitational_parameter,
                     double const planetarium_rotation_in_degrees);

// Deletes |plugin|.
extern "C" DLLEXPORT
void DestroyPlugin(Plugin* plugin);

// Calls |plugin->InsertCelestial| with the arguments given.
extern "C" DLLEXPORT
void InsertCelestial(Plugin* plugin, int const index,
                     double const gravitational_parameter,
                     int const parent,
                     XYZ const from_parent_position,
                     XYZ const from_parent_velocity);

// Calls |plugin->InsertCelestial| with the arguments given.
extern "C" DLLEXPORT
void UpdateCelestialHierarchy(Plugin* plugin, int const index,
                              int const parent);

// Calls |plugin->InsertOrKeepVessel| with the arguments given.
extern "C" DLLEXPORT
void InsertOrKeepVessel(Plugin* plugin, char const* guid, int const parent);

// Calls |plugin->SetVesselStateOffset| with the arguments given.
extern "C" DLLEXPORT
void SetVesselStateOffset(Plugin* plugin, char const* guid,
                          XYZ const from_parent_position,
                          XYZ const from_parent_velocity);

// Calls |plugin->VesselDisplacementFromParent| with the arguments given.
extern "C" DLLEXPORT
XYZ VesselDisplacementFromParent(Plugin* plugin,
                                 char const* guid);

// Calls |plugin->VesselParentRelativeVelocity| with the arguments given.
extern "C" DLLEXPORT
XYZ VesselParentRelativeVelocity(Plugin* plugin,
                                 char const* guid);

// Calls |plugin->CelestialDisplacementFromParent| with the arguments given.
extern "C" DLLEXPORT
XYZ CelestialDisplacementFromParent(Plugin* plugin, int const index);

// Calls |plugin->CelestialParentRelativeVelocity| with the arguments given. 
extern "C" DLLEXPORT
XYZ CelestialParentRelativeVelocity(Plugin* plugin, int const index);

extern "C" DLLEXPORT
char const* SayHello();

}  // namespace ksp_plugin
}  // namespace principia
