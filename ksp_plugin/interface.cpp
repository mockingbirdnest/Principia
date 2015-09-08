#include "ksp_plugin/interface.hpp"

#include <cctype>
#include <cstring>
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
using geometry::RadiusLatitudeLongitude;
using physics::MassiveBody;
using physics::OblateBody;
using quantities::Pow;
using si::AstronomicalUnit;
using si::Day;
using si::Degree;
using si::Kilo;
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

// Similar to std::stod, but uses LOG(FATAL) instead of exceptions.
double ParseDouble(std::string const& s, not_null<std::size_t*> size) {
  char* interpreted_end;
  char const* const c_string = s.c_str();
  double result = std::strtod(c_string, &interpreted_end);
  *size = interpreted_end - c_string;
  CHECK_GT(*size, 0) << "invalid floating-point number " << s;
  return result;
}

double ParseQuantity(std::string const& s, not_null<std::string*> unit) {
  std::size_t i;
  double magnitude = ParseDouble(s, &i);
  unit->clear();
  for (; i < s.length(); ++i) {
    if (!std::isspace(s[i])) {
      *unit += s[i];
    }
  }
  return magnitude;
}

Length ParseLength(std::string const& s) {
  std::string unit;
  double magnitude = ParseQuantity(s, &unit);
  if (unit == "m") {
    return magnitude * Metre;
  } else if (unit == "km") {
    return magnitude * Kilo(Metre);
  } else if (unit == "au") {
    return magnitude * AstronomicalUnit;
  } else {
    LOG(FATAL) << "unsupported unit of length " << unit;
    base::noreturn();
  }
}

Speed ParseSpeed(std::string const& s) {
  std::string unit;
  double magnitude = ParseQuantity(s, &unit);
  if (unit == "m/s") {
    return magnitude * Metre / Second;
  } else if (unit == "km/s") {
    return magnitude * Kilo(Metre) / Second;
  } else if (unit == "km/d") {
    return magnitude * Kilo(Metre) / Day;
  } else if (unit == "au/d") {
    return magnitude * AstronomicalUnit / Day;
  } else {
    LOG(FATAL) << "unsupported unit of speed " << unit;
    base::noreturn();
  }
}

Angle ParseAngle(std::string const& s) {
  std::string unit;
  double magnitude = ParseQuantity(s, &unit);
  if (unit == "deg" || unit == "°") {
    return magnitude * Degree;
  } else if (unit == "rad") {
    return magnitude * Radian;
  } else {
    LOG(FATAL) << "unsupported unit of angle " << unit;
    base::noreturn();
  }
}

GravitationalParameter ParseGravitationalParameter(std::string const& s) {
  std::string unit;
  double magnitude = ParseQuantity(s, &unit);
  if (unit == "m^3/s^2") {
    return magnitude * Pow<3>(Metre) / Pow<2>(Second);
  } else if (unit == "km^3/s^2") {
    return magnitude * Pow<3>(Kilo(Metre)) / Pow<2>(Second);
  } else if (unit == "km^3/d^2") {
    return magnitude * Pow<3>(Kilo(Metre)) / Pow<2>(Day);
  } else if (unit == "au^3/d^2") {
    return magnitude * Pow<3>(AstronomicalUnit) / Pow<2>(Day);
  } else {
    LOG(FATAL) << "unsupported unit of gravitational parameter " << unit;
    base::noreturn();
  }
}

double ParseDimensionless(std::string const& s) {
  std::string unit;
  double magnitude = ParseQuantity(s, &unit);
  CHECK(unit.empty()) << unit;
  return magnitude;
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

    google::protobuf::SetLogHandler(
        [](google::protobuf::LogLevel const level,
           char const* const filename,
           int const line,
           std::string const& message) {
          LOG_AT_LEVEL(level) << "[" << filename << ":" << line << "] "
                              << message;
        });

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
                             double const planetarium_rotation_in_degrees) {
  LOG(INFO) << "Constructing Principia plugin";
  not_null<std::unique_ptr<Plugin>> result = make_not_null_unique<Plugin>(
      Instant(initial_time * Second),
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

void principia__DirectlyInsertMassiveCelestial(
    Plugin* const plugin,
    int const celestial_index,
    int const* parent_index,
    char const* gravitational_parameter,
    char const* x,
    char const* y,
    char const* z,
    char const* vx,
    char const* vy,
    char const* vz) {
  CHECK_NOTNULL(plugin)->
      DirectlyInsertCelestial(
          celestial_index,
          parent_index,
          {Barycentric::origin +
               Displacement<Barycentric>({ParseLength(x),
                                          ParseLength(y),
                                          ParseLength(z)}),
               Velocity<Barycentric>({ParseSpeed(vx),
                                      ParseSpeed(vy),
                                      ParseSpeed(vz)})},
          std::make_unique<MassiveBody>(
              ParseGravitationalParameter(gravitational_parameter)));
}

void principia__DirectlyInsertOblateCelestial(
    Plugin* const plugin,
    int const celestial_index,
    int const* parent_index,
    char const* gravitational_parameter,
    char const* axis_right_ascension,
    char const* axis_declination,
    char const* j2,
    char const* reference_radius,
    char const* x,
    char const* y,
    char const* z,
    char const* vx,
    char const* vy,
    char const* vz) {
  CHECK_NOTNULL(plugin)->
    DirectlyInsertCelestial(
      celestial_index,
      parent_index,
      {Barycentric::origin +
       Displacement<Barycentric>({ParseLength(x),
                                  ParseLength(y),
                                  ParseLength(z)}),
       Velocity<Barycentric>({ParseSpeed(vx),
                              ParseSpeed(vy),
                              ParseSpeed(vz)})},
      std::make_unique<OblateBody<Barycentric>>(
          ParseGravitationalParameter(gravitational_parameter),
          ParseDimensionless(j2),
          ParseLength(reference_radius),
          Vector<double, Barycentric>(
              RadiusLatitudeLongitude(
                  1.0,
                  ParseAngle(axis_declination),
                  ParseAngle(axis_right_ascension)).ToCartesian())));
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

void principia__InsertSun(Plugin* const plugin,
                          int const celestial_index,
                          double const gravitational_parameter) {
  CHECK_NOTNULL(plugin)->InsertSun(
      celestial_index,
      gravitational_parameter * SIUnit<GravitationalParameter>());
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

int principia__ManœuvreCount(Plugin const* const plugin,
                             char const* const vessel_guid) {
  return CHECK_NOTNULL(plugin)->ManœuvreCount(vessel_guid);
}

Manœuvre<Barycentric> const* principia__VesselManœuvre(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  return &CHECK_NOTNULL(plugin)->VesselManœuvre(vessel_guid, index);
}

void principia__SetVesselManœuvre(Plugin const* const plugin,
                                  char const* const vessel_guid,
                                  int const index,
                                  Manœuvre<Barycentric> const** manœuvre) {
  CHECK_NOTNULL(plugin)->SetVesselManœuvre(vessel_guid,
                                           index,
                                           TakeOwnership(manœuvre));
}

void principia__InsertVesselManœuvre(Plugin const* const plugin,
                                     char const* const vessel_guid,
                                     int const index,
                                     Manœuvre<Barycentric> const** manœuvre) {
  CHECK_NOTNULL(plugin)->InsertVesselManœuvre(vessel_guid,
                                              index,
                                              TakeOwnership(manœuvre));
}

void principia__ClearVesselManœuvres(Plugin const* const plugin,
                                     char const* const vessel_guid) {
  CHECK_NOTNULL(plugin)->DeleteVesselManœuvres(vessel_guid);
}

XYZ principia__ManœuvreΔv(Plugin const* plugin,
                          Manœuvre<Barycentric> const* manœuvre) {
  return ToXYZ((CHECK_NOTNULL(plugin)->ManœuvreΔv(*CHECK_NOTNULL(manœuvre)) /
                    (Metre / Second)).coordinates());
}

void principia__UpdatePrediction(Plugin const* const plugin,
                                 char const* const vessel_guid) {
  CHECK_NOTNULL(plugin)->UpdatePrediction(vessel_guid);
}

void principia__UpdateFlightPlan(Plugin const* const plugin,
                                 char const* const vessel_guid,
                                 double const last_time) {
  CHECK_NOTNULL(plugin)->UpdateFlightPlan(vessel_guid,
                                          Instant(last_time * Second));
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

bool principia__HasPrediction(Plugin const* const plugin,
                              char const* const vessel_guid) {
  return CHECK_NOTNULL(plugin)->HasPrediction(vessel_guid);
}

LineAndIterator* principia__RenderedPrediction(
    Plugin* const plugin,
    char const* vessel_guid,
    RenderingTransforms* const transforms,
    XYZ const sun_world_position) {
  RenderedTrajectory<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedPrediction(
          vessel_guid,
          transforms,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return result.release();
}

int principia__FlightPlanSize(Plugin const* const plugin,
                              char const* const vessel_guid) {
  return CHECK_NOTNULL(plugin)->FlightPlanSize(vessel_guid);
}

LineAndIterator* principia__RenderedFlightPlan(
    Plugin* const plugin,
    char const* vessel_guid,
    int const plan_phase,
    RenderingTransforms* const transforms,
    XYZ const sun_world_position) {
  RenderedTrajectory<World> rendered_trajectory =
      CHECK_NOTNULL(plugin)->RenderedFlightPlan(
          vessel_guid,
          plan_phase,
          transforms,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return result.release();
}

void principia__set_prediction_length(Plugin* const plugin,
                                      double const t) {
  CHECK_NOTNULL(plugin)->set_prediction_length(t * Second);
}

void principia__set_prediction_length_tolerance(Plugin* const plugin,
                                                double const t) {
  CHECK_NOTNULL(plugin)->set_prediction_length_tolerance(t * Metre);
}

void principia__set_prediction_speed_tolerance(Plugin* const plugin,
                                               double const t) {
  CHECK_NOTNULL(plugin)->set_prediction_speed_tolerance(t * Metre / Second);
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
