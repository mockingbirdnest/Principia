#include "ksp_plugin/interface.hpp"

#include <string>

using principia::si::Degree;
using principia::si::Metre;

namespace principia {
namespace ksp_plugin {

void InitGoogleLogging() {
#ifdef _MSC_VER
  FILE* file;
  freopen_s(&file, "stderr.log", "w", stderr);
#else
  std::freopen("stderr.log", "w", stderr);
#endif
  google::SetStderrLogging(google::INFO);
  google::SetLogDestination(google::FATAL, "glog/Principia/FATAL.");
  google::SetLogDestination(google::ERROR, "glog/Principia/ERROR.");
  google::SetLogDestination(google::WARNING, "glog/Principia/WARNING.");
  google::SetLogDestination(google::INFO, "glog/Principia/INFO.");
  FLAGS_v = 1;
  // Buffer severities <= |INFO|, i.e., don't buffer.
  FLAGS_logbuflevel = google::INFO - 1;
  google::InitGoogleLogging("Principia");
  LOG(INFO) << "Initialized Google logging for Principia";
}

void LogInfo(char const* message) {
  LOG(INFO) << message;
}

void LogWarning(char const* message) {
  LOG(WARNING) << message;
}

void LogError(char const* message) {
  LOG(ERROR) << message;
}

void LogFatal(char const* message) {
  LOG(FATAL) << message;
}

Plugin* NewPlugin(double const initial_time,
                  int const sun_index,
                  double const sun_gravitational_parameter,
                  double const planetarium_rotation_in_degrees) {
  LOG(INFO) << "Constructing Principia plugin";
  return new Plugin(
      Instant(initial_time * Second),
      sun_index,
      sun_gravitational_parameter * SIUnit<GravitationalParameter>(),
      planetarium_rotation_in_degrees * Degree);
  LOG(INFO) << "Plugin constructed";
}

void DeletePlugin(Plugin const** const plugin) {
  LOG(INFO) << "Destroying Principia plugin";
  CHECK_NOTNULL(plugin);
  delete *plugin;
  *plugin = nullptr;
  LOG(INFO) << "Plugin destroyed";
}

// NOTE(egg): The |* (Metre / Second)| might be slower than |* SIUnit<Speed>()|,
// but it is more readable. This will be resolved once we have constexpr.

void InsertCelestial(Plugin* const plugin,
                     int const celestial_index,
                     double const gravitational_parameter,
                     int const parent_index,
                     XYZ const from_parent_position,
                     XYZ const from_parent_velocity) {
  CHECK_NOTNULL(plugin)->InsertCelestial(
      celestial_index,
      gravitational_parameter * SIUnit<GravitationalParameter>(),
      parent_index,
      Displacement<AliceSun>({from_parent_position.x * Metre,
                              from_parent_position.y * Metre,
                              from_parent_position.z * Metre}),
      Velocity<AliceSun>({from_parent_velocity.x * (Metre / Second),
                          from_parent_velocity.y * (Metre / Second),
                          from_parent_velocity.z * (Metre / Second)}));
}

void UpdateCelestialHierarchy(Plugin const* const plugin,
                              int const celestial_index,
                              int const parent_index) {
  CHECK_NOTNULL(plugin)->UpdateCelestialHierarchy(celestial_index,
                                                  parent_index);
}

void EndInitialization(Plugin* const plugin) {
  CHECK_NOTNULL(plugin)->EndInitialization();
}

void InsertOrKeepVessel(Plugin* const plugin,
                        char const* vessel_guid,
                        int const parent_index) {
  CHECK_NOTNULL(plugin)->InsertOrKeepVessel(vessel_guid, parent_index);
}

void SetVesselStateOffset(Plugin* const plugin,
                          char const* vessel_guid,
                          XYZ const from_parent_position,
                          XYZ const from_parent_velocity) {
  CHECK_NOTNULL(plugin)->SetVesselStateOffset(
      vessel_guid,
      Displacement<AliceSun>({from_parent_position.x * Metre,
                              from_parent_position.y * Metre,
                              from_parent_position.z * Metre}),
      Velocity<AliceSun>({from_parent_velocity.x * (Metre / Second),
                          from_parent_velocity.y * (Metre / Second),
                          from_parent_velocity.z * (Metre / Second)}));
}

void AdvanceTime(Plugin* const plugin,
                 double const t,
                 double const planetarium_rotation) {
  CHECK_NOTNULL(plugin)->AdvanceTime(Instant(t * Second),
                                     planetarium_rotation * Degree);
}

XYZ VesselDisplacementFromParent(Plugin const* const plugin,
                                 char const* vessel_guid) {
  R3Element<Length> const result =
      CHECK_NOTNULL(plugin)->
          VesselDisplacementFromParent(vessel_guid).coordinates();
  return {result.x / Metre, result.y / Metre, result.z / Metre};
}

XYZ VesselParentRelativeVelocity(Plugin const* const plugin,
                                 char const* vessel_guid) {
  R3Element<Speed> const result =
      CHECK_NOTNULL(plugin)->
          VesselParentRelativeVelocity(vessel_guid).coordinates();
  return {result.x / (Metre / Second),
          result.y / (Metre / Second),
          result.z / (Metre / Second)};
}

XYZ CelestialDisplacementFromParent(Plugin const* const plugin,
                                    int const celestial_index) {
  R3Element<Length> const result =
      CHECK_NOTNULL(plugin)->
          CelestialDisplacementFromParent(celestial_index).coordinates();
  return {result.x / Metre, result.y / Metre, result.z / Metre};;
}

XYZ CelestialParentRelativeVelocity(Plugin const* const plugin,
                                    int const celestial_index) {
  R3Element<Speed> const result =
      CHECK_NOTNULL(plugin)->
          CelestialParentRelativeVelocity(celestial_index).coordinates();
  return {result.x / (Metre / Second),
          result.y / (Metre / Second),
          result.z / (Metre / Second)};
}

char const* SayHello() {
  return "Hello from native C++!";
}

}  // namespace ksp_plugin
}  // namespace principia
