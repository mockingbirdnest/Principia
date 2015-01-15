#include "ksp_plugin/interface.hpp"

#include <string>
#include <utility>
#include <vector>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/version.hpp"
#include "ksp_plugin/part.hpp"

using principia::base::make_not_null_unique;
using principia::geometry::Displacement;
using principia::ksp_plugin::AliceSun;
using principia::ksp_plugin::LineSegment;
using principia::ksp_plugin::Part;
using principia::ksp_plugin::PartId;
using principia::ksp_plugin::RenderedTrajectory;
using principia::ksp_plugin::World;
using principia::quantities::Pow;
using principia::si::Degree;
using principia::si::Metre;
using principia::si::Second;
using principia::si::Tonne;

namespace {

// Takes ownership of |**pointer| and returns it to the caller.  Nulls
// |*pointer|.  |pointer| must not be null.  No transfer of ownership of
// |*pointer|.
template<typename T>
std::unique_ptr<T> TakeOwnership(T** const pointer) {
  CHECK_NOTNULL(pointer);
  std::unique_ptr<T> owned_pointer(*pointer);
  *pointer = nullptr;
  return owned_pointer;
}

R3Element<double> ToR3Element(XYZ const& xyz) {
  return {xyz.x, xyz.y, xyz.z};
}

XYZ ToXYZ(R3Element<double> const& r3_element) {
  return {r3_element.x, r3_element.y, r3_element.z};
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
    google::InitGoogleLogging("Principia");
    LOG(INFO) << "Initialized Google logging for Principia";
    LOG(INFO) << "Principia version " << principia::base::kVersion
              << " built on " << principia::base::kBuildDate
              << " by " << principia::base::kCompilerName
              << " version " << principia::base::kCompilerVersion
              << " for " << principia::base::kOperatingSystem
              << " " << principia::base::kArchitecture;
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

Plugin* principia__NewPlugin(double const initial_time,
                             int const sun_index,
                             double const sun_gravitational_parameter,
                             double const planetarium_rotation_in_degrees) {
  LOG(INFO) << "Constructing Principia plugin";
  not_null<std::unique_ptr<Plugin>> result = make_not_null_unique<Plugin>(
      Instant(initial_time * Second),
      sun_index,
      sun_gravitational_parameter * SIUnit<GravitationalParameter>(),
      planetarium_rotation_in_degrees * Degree);
  LOG(INFO) << "Plugin constructed";
  return result.release();
}

void principia__DeletePlugin(Plugin const** const plugin) {
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

void principia__InsertCelestial(Plugin* const plugin,
                                int const celestial_index,
                                double const gravitational_parameter,
                                int const parent_index,
                                QP const from_parent) {
  CHECK_NOTNULL(plugin)->InsertCelestial(
      celestial_index,
      gravitational_parameter * SIUnit<GravitationalParameter>(),
      parent_index,
      RelativeDegreesOfFreedom<AliceSun>(
          Displacement<AliceSun>(ToR3Element(from_parent.q) * Metre),
          Velocity<AliceSun>(ToR3Element(from_parent.p) * (Metre / Second))));
}

void principia__UpdateCelestialHierarchy(Plugin const* const plugin,
                                         int const celestial_index,
                                         int const parent_index) {
  CHECK_NOTNULL(plugin)->UpdateCelestialHierarchy(celestial_index,
                                                  parent_index);
}

void principia__EndInitialization(Plugin* const plugin) {
  CHECK_NOTNULL(plugin)->EndInitialization();
}

bool principia__InsertOrKeepVessel(Plugin* const plugin,
                                   char const* vessel_guid,
                                   int const parent_index) {
  return CHECK_NOTNULL(plugin)->InsertOrKeepVessel(vessel_guid, parent_index);
}

void principia__SetVesselStateOffset(Plugin* const plugin,
                                     char const* vessel_guid,
                                     QP const from_parent) {
  CHECK_NOTNULL(plugin)->SetVesselStateOffset(
      vessel_guid,
      RelativeDegreesOfFreedom<AliceSun>(
          Displacement<AliceSun>(ToR3Element(from_parent.q) * Metre),
          Velocity<AliceSun>(ToR3Element(from_parent.p) * (Metre / Second))));
}

void principia__AdvanceTime(Plugin* const plugin,
                            double const t,
                            double const planetarium_rotation) {
  CHECK_NOTNULL(plugin)->AdvanceTime(Instant(t * Second),
                                     planetarium_rotation * Degree);
}

QP principia__VesselFromParent(Plugin const* const plugin,
                               char const* vessel_guid) {
  RelativeDegreesOfFreedom<AliceSun> const result =
      CHECK_NOTNULL(plugin)->VesselFromParent(vessel_guid);
  return {ToXYZ(result.displacement().coordinates() / Metre),
          ToXYZ(result.velocity().coordinates() / (Metre / Second))};
}

QP principia__CelestialFromParent(Plugin const* const plugin,
                                   int const celestial_index) {
  RelativeDegreesOfFreedom<AliceSun> const result =
      CHECK_NOTNULL(plugin)->CelestialFromParent(celestial_index);
  return {ToXYZ(result.displacement().coordinates() / Metre),
          ToXYZ(result.velocity().coordinates() / (Metre / Second))};
}

Transforms<Barycentric, Rendering, Barycentric>*
principia__NewBodyCentredNonRotatingTransforms(Plugin const* const plugin,
                                               int const reference_body_index) {
  return CHECK_NOTNULL(plugin)->
      NewBodyCentredNonRotatingTransforms(reference_body_index).release();
}

Transforms<Barycentric, Rendering, Barycentric>*
principia__NewBarycentricRotatingTransforms(Plugin const* const plugin,
                                            int const primary_index,
                                            int const secondary_index) {
  return CHECK_NOTNULL(plugin)->
      NewBarycentricRotatingTransforms(
          primary_index, secondary_index).release();
}

void principia__DeleteTransforms(
    Transforms<Barycentric, Rendering, Barycentric>** const transforms) {
  TakeOwnership(transforms);
}

LineAndIterator* principia__RenderedVesselTrajectory(
    Plugin const* const plugin,
    char const* vessel_guid,
    Transforms<Barycentric, Rendering, Barycentric>* const transforms,
    XYZ const sun_world_position) {
  RenderedTrajectory<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedVesselTrajectory(
          vessel_guid,
          transforms,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return result.release();
}

int principia__NumberOfSegments(LineAndIterator const* line_and_iterator) {
  return CHECK_NOTNULL(line_and_iterator)->rendered_trajectory.size();
}

XYZSegment principia__FetchAndIncrement(
    LineAndIterator* const line_and_iterator) {
  CHECK_NOTNULL(line_and_iterator);
  CHECK(line_and_iterator->it != line_and_iterator->rendered_trajectory.end());
  LineSegment<World> const result = *line_and_iterator->it;
  ++line_and_iterator->it;
  return {ToXYZ((result.begin - World::origin).coordinates() / Metre),
          ToXYZ((result.end - World::origin).coordinates() / Metre)};
}

bool principia__AtEnd(LineAndIterator* const line_and_iterator) {
  CHECK_NOTNULL(line_and_iterator);
  return line_and_iterator->it == line_and_iterator->rendered_trajectory.end();
}

void principia__DeleteLineAndIterator(
    LineAndIterator** const line_and_iterator) {
  TakeOwnership(line_and_iterator);
}

XYZ principia__VesselWorldPosition(Plugin const* const plugin,
                                   char const* vessel_guid,
                                   XYZ const parent_world_position) {
  Position<World> const result = CHECK_NOTNULL(plugin)->VesselWorldPosition(
      vessel_guid,
      World::origin + Displacement<World>(
                          ToR3Element(parent_world_position) * Metre));
  return ToXYZ((result - World::origin).coordinates() / Metre);
}

XYZ principia__VesselWorldVelocity(Plugin const* const plugin,
                                   char const* vessel_guid,
                                   XYZ const parent_world_velocity,
                                   double const parent_rotation_period) {
  Velocity<World> const result = CHECK_NOTNULL(plugin)->VesselWorldVelocity(
      vessel_guid,
      Velocity<World>(ToR3Element(parent_world_velocity) * (Metre / Second)),
      parent_rotation_period * Second);
  return ToXYZ(result.coordinates() / (Metre / Second));
}

void principia__AddVesselToNextPhysicsBubble(Plugin* const plugin,
                                             char const* vessel_guid,
                                             KSPPart const* const parts,
                                             int count) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(count);
  std::vector<principia::ksp_plugin::IdAndOwnedPart> vessel_parts;
  vessel_parts.reserve(count);
  for (KSPPart const* part = parts; part < parts + count; ++part) {
    vessel_parts.push_back(
        std::make_pair(
            part->id,
            make_not_null_unique<Part<World>>(
                DegreesOfFreedom<World>(
                    World::origin +
                        Displacement<World>(
                            ToR3Element(part->world_position) * Metre),
                    Velocity<World>(
                        ToR3Element(part->world_velocity) * (Metre / Second))),
                part->mass * Tonne,
                Vector<Acceleration, World>(
                    ToR3Element(
                        part->gravitational_acceleration_to_be_applied_by_ksp) *
                    (Metre / Pow<2>(Second))))));
  }
  CHECK_NOTNULL(plugin)->AddVesselToNextPhysicsBubble(vessel_guid,
                                                      std::move(vessel_parts));
}

XYZ principia__BubbleDisplacementCorrection(Plugin const* const plugin,
                                            XYZ const sun_position) {
  Displacement<World> const result =
      CHECK_NOTNULL(plugin)->BubbleDisplacementCorrection(
          World::origin + Displacement<World>(
                              ToR3Element(sun_position) * Metre));
  return ToXYZ(result.coordinates() / Metre);
}

XYZ principia__BubbleVelocityCorrection(Plugin const* const plugin,
                                        int const reference_body_index) {
  Velocity<World> const result =
      CHECK_NOTNULL(plugin)->BubbleVelocityCorrection(reference_body_index);
  return ToXYZ(result.coordinates() / (Metre / Second));
}

double principia__current_time(Plugin const* const plugin) {
  return (CHECK_NOTNULL(plugin)->current_time() - Instant()) / Second;
}

char const* principia__SayHello() {
  return "Hello from native C++!";
}
