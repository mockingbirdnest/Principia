#pragma once

#include <type_traits>

#include "ksp_plugin/plugin.hpp"

// DLL-exported functions for interfacing with Platform Invocation Services.

#if defined(DLLEXPORT)
#error "DLLEXPORT already defined"
#else
#if defined(_WIN32)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __attribute__((visibility("default")))
#endif
#endif

namespace principia {
namespace ksp_plugin {

extern "C"
struct XYZ {
  double x, y, z;
};

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
// This will always evaluate its argument even if the corresponding log severity
// is disabled, so it is less efficient than LOG(INFO). It will not report the
// line and file of the caller.
extern "C" DLLEXPORT
void LogInfo(char const* message);
extern "C" DLLEXPORT
void LogWarning(char const* message);
extern "C" DLLEXPORT
void LogError(char const* message);
extern "C" DLLEXPORT
void LogFatal(char const* message);

// Returns a pointer to a plugin constructed with the arguments given.
// The caller takes ownership of the result.
extern "C" DLLEXPORT
Plugin* NewPlugin(double const initial_time,
                  int const sun_index,
                  double const sun_gravitational_parameter,
                  double const planetarium_rotation_in_degrees);

// Deletes and nulls |*plugin|.
// No transfer of ownership of |*plugin|, takes ownership of |**plugin|.
extern "C" DLLEXPORT
void DeletePlugin(Plugin const** const plugin);

// Calls |plugin->InsertCelestial| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
void InsertCelestial(Plugin* const plugin,
                     int const celestial_index,
                     double const gravitational_parameter,
                     int const parent_index,
                     XYZ const from_parent_position,
                     XYZ const from_parent_velocity);

// Calls |plugin->UpdateCelestialHierarchy| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
void UpdateCelestialHierarchy(Plugin const* const plugin,
                              int const celestial_index,
                              int const parent_index);

// Calls |plugin->InsertOrKeepVessel| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
void InsertOrKeepVessel(Plugin* const plugin,
                        char const* vessel_guid,
                        int const parent_index);

// Calls |plugin->SetVesselStateOffset| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
void SetVesselStateOffset(Plugin const* const plugin,
                          char const* vessel_guid,
                          XYZ const from_parent_position,
                          XYZ const from_parent_velocity);

// Calls |plugin->VesselDisplacementFromParent| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
XYZ VesselDisplacementFromParent(Plugin const* const plugin,
                                 char const* vessel_guid);

// Calls |plugin->VesselParentRelativeVelocity| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
XYZ VesselParentRelativeVelocity(Plugin const* const plugin,
                                 char const* vessel_guid);

// Calls |plugin->CelestialDisplacementFromParent| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
XYZ CelestialDisplacementFromParent(Plugin const* const plugin,
                                    int const celestial_index);

// Calls |plugin->CelestialParentRelativeVelocity| with the arguments given.
// |plugin| should not be null. No transfer of ownership.
extern "C" DLLEXPORT
XYZ CelestialParentRelativeVelocity(Plugin const* const plugin,
                                    int const celestial_index);

// Says hello, convenient for checking that calls to the dll work.
extern "C" DLLEXPORT
char const* SayHello();

}  // namespace ksp_plugin
}  // namespace principia

#undef DLLEXPORT
