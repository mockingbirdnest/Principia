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
  google::InitGoogleLogging("Principia");
  LOG(INFO) << "Initialized Google logging for Principia.";
}

void LOGINFO(char const* message) {
  LOG(INFO) << message;
}

void LOGWARNING(char const* message) {
  LOG(WARNING) << message;
}

void LOGERROR(char const* message) {
  LOG(ERROR) << message;
}

void LOGFATAL(char const* message) {
  LOG(FATAL) << message;
}

Plugin* CreatePlugin(double const initial_time, int const sun_index,
                     double const sun_gravitational_parameter,
                     double const planetarium_rotation_in_degrees) {
  LOG(INFO) << "Creating Principia...";
  return new Plugin(
      Instant(initial_time * Second),
      sun_index,
      sun_gravitational_parameter * SIUnit<GravitationalParameter>(),
      planetarium_rotation_in_degrees * Degree);
  LOG(INFO) << "Plugin created.";
}

void DestroyPlugin(Plugin* plugin) {
  delete plugin;
  plugin = nullptr;
  LOG(INFO) << "Destroyed Principia.";
}

void InsertCelestial(Plugin* plugin, int const index,
                     double const gravitational_parameter,
                     int const parent,
                     XYZ const from_parent_position,
                     XYZ const from_parent_velocity) {
  plugin->InsertCelestial(
      index,
      gravitational_parameter * SIUnit<GravitationalParameter>(),
      parent,
      Displacement<AliceSun>({from_parent_position.x * Metre,
                              from_parent_position.y * Metre,
                              from_parent_position.z * Metre}),
      Velocity<AliceSun>({from_parent_velocity.x * Metre / Second,
                          from_parent_velocity.y * Metre / Second,
                          from_parent_velocity.z * Metre / Second}));
}

void UpdateCelestialHierarchy(Plugin* plugin, int const index,
                              int const parent) {
  plugin->UpdateCelestialHierarchy(index, parent);
}

void InsertOrKeepVessel(Plugin* plugin, char const* guid, int const parent) {
  plugin->InsertOrKeepVessel(guid, parent);
}

void SetVesselStateOffset(Plugin* plugin, char const* guid,
                          XYZ const from_parent_position,
                          XYZ const from_parent_velocity) {
  plugin->SetVesselStateOffset(
      guid,
      Displacement<AliceSun>({from_parent_position.x * Metre,
                              from_parent_position.y * Metre,
                              from_parent_position.z * Metre}),
      Velocity<AliceSun>({from_parent_velocity.x * Metre / Second,
                          from_parent_velocity.y * Metre / Second,
                          from_parent_velocity.z * Metre / Second}));
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
