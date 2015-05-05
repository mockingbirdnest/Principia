#include "ksp_plugin/interface.hpp"

#include <string>
#include <utility>
#include <vector>
#if OS_WIN
#define NOGDI
#include <windows.h>
#include <psapi.h>
#endif

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "base/version.hpp"
#include "ksp_plugin/part.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using base::Bytes;
using base::HexadecimalDecode;
using base::HexadecimalEncode;
using base::make_not_null_unique;
using base::PullSerializer;
using base::PushDeserializer;
using base::UniqueBytes;
using geometry::Displacement;
using geometry::Quaternion;
using quantities::Pow;
using si::Degree;
using si::Metre;
using si::Second;
using si::Tonne;

namespace ksp_plugin {

namespace {

int const kChunkSize = 64 << 10;
int const kNumberOfChunks = 8;

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

template<typename T>
std::unique_ptr<T[]> TakeOwnershipArray(T** const pointer) {
  CHECK_NOTNULL(pointer);
  std::unique_ptr<T[]> owned_pointer(*pointer);
  *pointer = nullptr;
  return owned_pointer;
}

R3Element<double> ToR3Element(XYZ const& xyz) {
  return {xyz.x, xyz.y, xyz.z};
}

XYZ ToXYZ(R3Element<double> const& r3_element) {
  return {r3_element.x, r3_element.y, r3_element.z};
}

WXYZ ToWXYZ(Quaternion const& quaternion) {
  return {quaternion.real_part(),
          quaternion.imaginary_part().x,
          quaternion.imaginary_part().y,
          quaternion.imaginary_part().z};
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
    google::SetLogDestination(google::FATAL, "glog/Principia/FATAL.");
    google::SetLogDestination(google::ERROR, "glog/Principia/ERROR.");
    google::SetLogDestination(google::WARNING, "glog/Principia/WARNING.");
    google::SetLogDestination(google::INFO, "glog/Principia/INFO.");
    google::InitGoogleLogging("Principia");
    LOG(INFO) << "Initialized Google logging for Principia";
    LOG(INFO) << "Principia version " << principia::base::kVersion
              << " built on " << principia::base::kBuildDate
              << " by " << principia::base::kCompilerName
              << " version " << principia::base::kCompilerVersion
              << " for " << principia::base::kOperatingSystem
              << " " << principia::base::kArchitecture;
#if OS_WIN
  MODULEINFO module_info;
  memset(&module_info, 0, sizeof(module_info));
  CHECK(GetModuleInformation(GetCurrentProcess(),
                             GetModuleHandle(TEXT("principia")),
                             &module_info,
                             sizeof(module_info)));
  LOG(INFO) << "Base address is " << module_info.lpBaseOfDll;
#endif
  }
}

void principia__SetBufferedLogging(int const max_severity) {
  FLAGS_logbuflevel = max_severity;
}

int principia__GetBufferedLogging() {
  return FLAGS_logbuflevel;
}

void principia__SetBufferDuration(int const seconds) {
  FLAGS_logbufsecs = seconds;
}

int principia__GetBufferDuration() {
  return FLAGS_logbufsecs;
}

void principia__SetSuppressedLogging(int const min_severity) {
  FLAGS_minloglevel = min_severity;
}

int principia__GetSuppressedLogging() {
  return FLAGS_minloglevel;
}

void principia__SetVerboseLogging(int const level) {
  FLAGS_v = level;
}

int principia__GetVerboseLogging() {
  return FLAGS_v;
}

void principia__SetStderrLogging(int const min_severity) {
  // NOTE(egg): We could use |FLAGS_stderrthreshold| instead, the difference
  // seems to be a mutex.
  google::SetStderrLogging(min_severity);
}

int principia__GetStderrLogging() {
  return FLAGS_stderrthreshold;
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

void principia__ForgetAllHistoriesBefore(Plugin* const plugin,
                                         double const t) {
  CHECK_NOTNULL(plugin)->ForgetAllHistoriesBefore(Instant(t * Second));
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

RenderingTransforms* principia__NewBodyCentredNonRotatingTransforms(
    Plugin const* const plugin,
    int const reference_body_index) {
  return CHECK_NOTNULL(plugin)->
      NewBodyCentredNonRotatingTransforms(reference_body_index).release();
}

RenderingTransforms* principia__NewBarycentricRotatingTransforms(
    Plugin const* const plugin,
    int const primary_index,
    int const secondary_index) {
  return CHECK_NOTNULL(plugin)->
      NewBarycentricRotatingTransforms(
          primary_index, secondary_index).release();
}

void principia__DeleteTransforms(RenderingTransforms** const transforms) {
  TakeOwnership(transforms);
}

LineAndIterator* principia__RenderedVesselTrajectory(
    Plugin const* const plugin,
    char const* vessel_guid,
    RenderingTransforms* const transforms,
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

LineAndIterator* principia__RenderedPrediction(
    Plugin* const plugin,
    RenderingTransforms* const transforms,
    XYZ const sun_world_position) {
  RenderedTrajectory<World> rendered_trajectory =
      CHECK_NOTNULL(plugin)->RenderedPrediction(
          transforms,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return result.release();
}

void principia__set_predicted_vessel(Plugin* const plugin,
                                     char const* vessel_guid) {
  CHECK_NOTNULL(plugin)->set_predicted_vessel(vessel_guid);
}

void principia__clear_predicted_vessel(Plugin* const plugin) {
  CHECK_NOTNULL(plugin)->clear_predicted_vessel();
}

void principia__set_prediction_length(Plugin* const plugin,
                                      double const t) {
  CHECK_NOTNULL(plugin)->set_prediction_length(t * Second);
}

void principia__set_prediction_step(Plugin* const plugin,
                                    double const t) {
  CHECK_NOTNULL(plugin)->set_prediction_step(t * Second);
}

bool principia__has_vessel(Plugin* const plugin,
                           char const* vessel_guid) {
  return CHECK_NOTNULL(plugin)->has_vessel(vessel_guid);
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

bool principia__PhysicsBubbleIsEmpty(Plugin const* const plugin) {
  return CHECK_NOTNULL(plugin)->PhysicsBubbleIsEmpty();
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

WXYZ principia__NavballOrientation(
    Plugin const* const plugin,
    RenderingTransforms* const transforms,
    XYZ const sun_world_position,
    XYZ const ship_world_position) {
  FrameField<World> const frame_field = CHECK_NOTNULL(plugin)->Navball(
      transforms,
      World::origin +
          Displacement<World>(ToR3Element(sun_world_position) * Metre));
  return ToWXYZ(
      frame_field(
          World::origin +
              Displacement<World>(
                  ToR3Element(ship_world_position) * Metre)).quaternion());
}

XYZ principia__VesselTangent(Plugin const* const plugin,
                             char const* vessel_guid,
                             RenderingTransforms* const transforms) {
  return ToXYZ(CHECK_NOTNULL(plugin)->
                   VesselTangent(vessel_guid, transforms).coordinates());
}

double principia__current_time(Plugin const* const plugin) {
  return (CHECK_NOTNULL(plugin)->current_time() - Instant()) / Second;
}

char const* principia__SerializePlugin(Plugin const* const plugin,
                                       PullSerializer** const serializer) {
  LOG(INFO) << __FUNCTION__;
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(serializer);

  // Create and start a serializer if the caller didn't provide one.
  if (*serializer == nullptr) {
    *serializer = new PullSerializer(kChunkSize, kNumberOfChunks);
    auto message = make_not_null_unique<serialization::Plugin>();
    plugin->WriteToMessage(message.get());
    (*serializer)->Start(std::move(message));
  }

  // Pull a chunk.
  Bytes bytes;
  bytes = (*serializer)->Pull();

  // If this is the end of the serialization, delete the serializer and return a
  // nullptr.
  if (bytes.size == 0) {
    TakeOwnership(serializer);
    return nullptr;
  }

  // Convert to hexadecimal and return to the client.
  std::int64_t const hexadecimal_size = (bytes.size << 1) + 1;
  UniqueBytes hexadecimal(hexadecimal_size);
  HexadecimalEncode(bytes, hexadecimal.get());
  hexadecimal.data.get()[hexadecimal_size - 1] = '\0';
  return reinterpret_cast<char const*>(hexadecimal.data.release());
}

void principia__DeletePluginSerialization(char const** const serialization) {
  LOG(INFO) << __FUNCTION__;
  TakeOwnershipArray(reinterpret_cast<uint8_t const**>(serialization));
}

void principia__DeserializePlugin(char const* const serialization,
                                  int const serialization_size,
                                  PushDeserializer** const deserializer,
                                  Plugin const** const plugin) {
  LOG(INFO) << __FUNCTION__;
  CHECK_NOTNULL(serialization);
  CHECK_NOTNULL(deserializer);
  CHECK_NOTNULL(plugin);

  // Create and start a deserializer if the caller didn't provide one.
  if (*deserializer == nullptr) {
    *deserializer = new PushDeserializer(kChunkSize, kNumberOfChunks);
    auto message = make_not_null_unique<serialization::Plugin>();
    (*deserializer)->Start(
        std::move(message),
        [plugin](google::protobuf::Message const& message) {
          *plugin = Plugin::ReadFromMessage(
              static_cast<serialization::Plugin const&>(message)).release();
        });
  }

  // Decode the hexadecimal representation.
  uint8_t const* const hexadecimal =
      reinterpret_cast<uint8_t const*>(serialization);
  int const hexadecimal_size = serialization_size;
  int const byte_size = hexadecimal_size >> 1;
  // Ownership of the following pointer is transfered to the deserializer using
  // the callback to |Push|.
  std::uint8_t* bytes = new uint8_t[byte_size];
  HexadecimalDecode({hexadecimal, hexadecimal_size}, {bytes, byte_size});

  // Push the data, taking ownership of it.
  (*deserializer)->Push(Bytes(&bytes[0], byte_size),
                        [bytes]() { delete bytes; });

  // If the data was empty, delete the deserializer.  This ensures that
  // |*plugin| is filled.
  if (byte_size == 0) {
    delete *deserializer;
  }
}

char const* principia__SayHello() {
  return "Hello from native C++!";
}

}  // namespace ksp_plugin
}  // namespace principia
