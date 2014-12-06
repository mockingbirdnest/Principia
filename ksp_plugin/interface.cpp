#include "ksp_plugin/interface.hpp"

#include <string>

#include "base/version.hpp"

using principia::geometry::Displacement;
using principia::ksp_plugin::AliceSun;
using principia::ksp_plugin::LineSegment;
using principia::ksp_plugin::RenderedTrajectory;
using principia::ksp_plugin::World;
using principia::si::Degree;
using principia::si::Metre;
using principia::si::Second;

namespace {

// Takes ownership of |**pointer| and returns it to the caller.  Nulls
// |*pointer|.  |pointer| must not be null.  No transfer of ownership of
// |*pointer|.
template<typename T>
std::unique_ptr<T> TakeOwnership(T** const pointer) {
  CHECK_NOTNULL(pointer);
  std::unique_ptr<T const> owned_pointer(*pointer);
  *pointer = nullptr;
  return owned_pointer;
}

}  // namespace

void principia__InitGoogleLogging() {
  if (google::IsGoogleLoggingInitialized()) {
    LOG(INFO) << "Google logging was already initialized, no action taken";
  } else {
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
    google::principia__InitGoogleLogging("Principia");
    LOG(INFO) << "Initialized Google logging for Principia";
    LOG(INFO) << "Principia version " << principia::base::kVersion
              << " built on " << principia::base::kBuildDate;
    // TODO(egg): by (compiler) for (ARCH, OS).
  }
}

void principia__LogInfo(char const* message) {
  LOG(INFO) << message;
}

void principia__LogWarning(char const* message) {
  LOG(WARNING) << message;
}

void principia__LogError(char const* message) {
  LOG(ERROR) << message;
}

void principia__LogFatal(char const* message) {
  LOG(FATAL) << message;
}

Plugin* NewPlugin(double const initial_time,
                  int const sun_index,
                  double const sun_gravitational_parameter,
                  double const planetarium_rotation_in_degrees) {
  LOG(INFO) << "Constructing Principia plugin";
  std::unique_ptr<Plugin> result = std::make_unique<Plugin>(
      Instant(initial_time * Second),
      sun_index,
      sun_gravitational_parameter * SIUnit<GravitationalParameter>(),
      planetarium_rotation_in_degrees * Degree);
  LOG(INFO) << "Plugin constructed";
  return result.release();
}

void DeletePlugin(Plugin const** const plugin) {
  LOG(INFO) << "Destroying Principia plugin";
  // We want to log before and after destroying the plugin since it is a pretty
  // significant event, so we take ownership inside a block.
  {
    TakeOwnership(plugin);
  }
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

bool InsertOrKeepVessel(Plugin* const plugin,
                        char const* vessel_guid,
                        int const parent_index) {
  return CHECK_NOTNULL(plugin)->InsertOrKeepVessel(vessel_guid, parent_index);
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
  return {result.x / Metre, result.y / Metre, result.z / Metre};
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

BodyCentredNonRotatingFrame const* NewBodyCentredNonRotatingFrame(
    Plugin const* const plugin,
    int const reference_body_index) {
  return CHECK_NOTNULL(plugin)->
      NewBodyCentredNonRotatingFrame(reference_body_index).release();
}

BarycentricRotatingFrame const* NewBarycentricRotatingFrame(
    Plugin const* const plugin,
    int const primary_index,
    int const secondary_index) {
  return CHECK_NOTNULL(plugin)->
      NewBarycentricRotatingFrame(primary_index, secondary_index).release();
}

void DeleteRenderingFrame(RenderingFrame const** const frame) {
  TakeOwnership(frame);
}

LineAndIterator* RenderedVesselTrajectory(Plugin const* const plugin,
                                          char const* vessel_guid,
                                          RenderingFrame const* frame,
                                          XYZ const sun_world_position) {
  RenderedTrajectory<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedVesselTrajectory(
          vessel_guid,
          *CHECK_NOTNULL(frame),
          World::origin + Displacement<World>({sun_world_position.x * Metre,
                                               sun_world_position.y * Metre,
                                               sun_world_position.z * Metre}));
  std::unique_ptr<LineAndIterator> result =
      std::make_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return result.release();
}

int NumberOfSegments(LineAndIterator const* line_and_iterator) {
  return CHECK_NOTNULL(line_and_iterator)->rendered_trajectory.size();
}

XYZSegment FetchAndIncrement(LineAndIterator* const line_and_iterator) {
  CHECK_NOTNULL(line_and_iterator);
  CHECK(line_and_iterator->it != line_and_iterator->rendered_trajectory.end());
  LineSegment<World> const result = *line_and_iterator->it;
  ++line_and_iterator->it;
  R3Element<Length> const begin = (result.begin - World::origin).coordinates();
  R3Element<Length> const end = (result.end - World::origin).coordinates();
  return {XYZ{begin.x / Metre, begin.y / Metre, begin.z / Metre},
          XYZ{end.x / Metre, end.y / Metre, end.z / Metre}};
}

bool AtEnd(LineAndIterator* const line_and_iterator) {
  CHECK_NOTNULL(line_and_iterator);
  return line_and_iterator->it == line_and_iterator->rendered_trajectory.end();
}

void DeleteLineAndIterator(LineAndIterator const** const line_and_iterator) {
  TakeOwnership(line_and_iterator);
}

XYZ VesselWorldPosition(Plugin const* const plugin,
                        char const* vessel_guid,
                        XYZ const parent_world_position) {
  Position<World> result = CHECK_NOTNULL(plugin)->VesselWorldPosition(
      vessel_guid,
      World::origin + Displacement<World>({parent_world_position.x * Metre,
                                           parent_world_position.y * Metre,
                                           parent_world_position.z * Metre}));
  R3Element<Length> const& coordinates = (result - World::origin).coordinates();
  return XYZ{coordinates.x / Metre,
             coordinates.y / Metre,
             coordinates.z / Metre};
}

XYZ VesselWorldVelocity(Plugin const* const plugin,
                        char const* vessel_guid,
                        XYZ const parent_world_velocity,
                        double const parent_rotation_period) {
  Velocity<World> result = CHECK_NOTNULL(plugin)->VesselWorldVelocity(
      vessel_guid,
      Velocity<World>({parent_world_velocity.x * Metre / Second,
                       parent_world_velocity.y * Metre / Second,
                       parent_world_velocity.z * Metre / Second}),
      parent_rotation_period * Second);
  R3Element<Speed> const& coordinates = result.coordinates();
  return XYZ{coordinates.x / (Metre / Second),
             coordinates.y / (Metre / Second),
             coordinates.z / (Metre / Second)};
}

char const* SayHello() {
  return "Hello from native C++!";
}
