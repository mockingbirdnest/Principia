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
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "journal/recorder.hpp"
#include "ksp_plugin/part.hpp"
#include "physics/solar_system.hpp"
#include "quantities/parser.hpp"
#include "serialization/astronomy.pb.h"
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
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::ParseQuantity;
using quantities::Pow;
using quantities::si::AstronomicalUnit;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using quantities::si::Tonne;

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
  journal::Method<journal::InitGoogleLogging> m;
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
  return m.Return();
}

void principia__ActivateRecorder(bool const activate) {
  // NOTE: Do not journal!  You'd end up with half a message in the journal and
  // that would cause trouble.
  if (activate && !journal::Recorder::IsActivated()) {
    // Build a name somewhat similar to that of the log files.
    auto const now = std::chrono::system_clock::now();
    std::time_t const time = std::chrono::system_clock::to_time_t(now);
    std::tm* const localtime = std::localtime(&time);
    std::stringstream name;
    name << std::put_time(localtime, "JOURNAL.%Y%m%d-%H%M%S");
    journal::Recorder* const recorder =
        new journal::Recorder(std::experimental::filesystem::path("glog") /
                              "Principia" / name.str());
    journal::Recorder::Activate(recorder);
  } else if (!activate && journal::Recorder::IsActivated()) {
    journal::Recorder::Deactivate();
  }
}

void principia__SetBufferedLogging(int const max_severity) {
  journal::Method<journal::SetBufferedLogging> m({max_severity});
  FLAGS_logbuflevel = max_severity;
  return m.Return();
}

int principia__GetBufferedLogging() {
  journal::Method<journal::GetBufferedLogging> m;
  return m.Return(FLAGS_logbuflevel);
}

void principia__SetBufferDuration(int const seconds) {
  journal::Method<journal::SetBufferDuration> m({seconds});
  FLAGS_logbufsecs = seconds;
  return m.Return();
}

int principia__GetBufferDuration() {
  journal::Method<journal::GetBufferDuration> m;
  return m.Return(FLAGS_logbufsecs);
}

void principia__SetSuppressedLogging(int const min_severity) {
  journal::Method<journal::SetSuppressedLogging> m({min_severity});
  FLAGS_minloglevel = min_severity;
  return m.Return();
}

int principia__GetSuppressedLogging() {
  journal::Method<journal::GetSuppressedLogging> m;
  return m.Return(FLAGS_minloglevel);
}

void principia__SetVerboseLogging(int const level) {
  journal::Method<journal::SetVerboseLogging> m({level});
  FLAGS_v = level;
  return m.Return();
}

int principia__GetVerboseLogging() {
  journal::Method<journal::GetVerboseLogging> m;
  return m.Return(FLAGS_v);
}

void principia__SetStderrLogging(int const min_severity) {
  journal::Method<journal::SetStderrLogging> m({min_severity});
  // NOTE(egg): We could use |FLAGS_stderrthreshold| instead, the difference
  // seems to be a mutex.
  google::SetStderrLogging(min_severity);
  return m.Return();
}

int principia__GetStderrLogging() {
  journal::Method<journal::GetStderrLogging> m;
  return m.Return(FLAGS_stderrthreshold);
}

void principia__LogInfo(char const* const text) {
  journal::Method<journal::LogInfo> m({text});
  LOG(INFO) << text;
  return m.Return();
}

void principia__LogWarning(char const* const text) {
  journal::Method<journal::LogWarning> m({text});
  LOG(WARNING) << text;
  return m.Return();
}

void principia__LogError(char const* const text) {
  journal::Method<journal::LogError> m({text});
  LOG(ERROR) << text;
  return m.Return();
}

void principia__LogFatal(char const* const text) {
  journal::Method<journal::LogFatal> m({text});
  LOG(FATAL) << text;
  return m.Return();
}

Plugin* principia__NewPlugin(double const initial_time,
                             double const planetarium_rotation_in_degrees) {
  journal::Method<journal::NewPlugin> m({initial_time,
                                         planetarium_rotation_in_degrees});
  LOG(INFO) << "Constructing Principia plugin";
  Instant const t0;
  not_null<std::unique_ptr<Plugin>> result = make_not_null_unique<Plugin>(
      t0 + initial_time * Second,
      planetarium_rotation_in_degrees * Degree);
  LOG(INFO) << "Plugin constructed";
  return m.Return(result.release());
}

void principia__DeletePlugin(Plugin const** const plugin) {
  CHECK_NOTNULL(plugin);
  journal::Method<journal::DeletePlugin> m({plugin}, {plugin});
  LOG(INFO) << "Destroying Principia plugin";
  // We want to log before and after destroying the plugin since it is a pretty
  // significant event, so we take ownership inside a block.
  {
    TakeOwnership(plugin);
  }
  LOG(INFO) << "Plugin destroyed";
  return m.Return();
}

void principia__DirectlyInsertCelestial(
    Plugin* const plugin,
    int const celestial_index,
    int const* const parent_index,
    char const* const gravitational_parameter,
    char const* const axis_right_ascension,
    char const* const axis_declination,
    char const* const j2,
    char const* const reference_radius,
    char const* const x,
    char const* const y,
    char const* const z,
    char const* const vx,
    char const* const vy,
    char const* const vz) {
  journal::Method<journal::DirectlyInsertCelestial> m({plugin,
                                                       celestial_index,
                                                       parent_index,
                                                       gravitational_parameter,
                                                       axis_right_ascension,
                                                       axis_declination,
                                                       j2,
                                                       reference_radius,
                                                       x, y, z,
                                                       vx, vy, vz});
  serialization::GravityModel::Body gravity_model;
  serialization::InitialState::Body initial_state;
  gravity_model.set_gravitational_parameter(gravitational_parameter);
  if (axis_right_ascension != nullptr) {
    gravity_model.set_axis_right_ascension(axis_right_ascension);
  }
  if (axis_declination != nullptr) {
    gravity_model.set_axis_declination(axis_declination);
  }
  if (j2 != nullptr) {
    gravity_model.set_j2(j2);
  }
  if (reference_radius != nullptr) {
    gravity_model.set_reference_radius(reference_radius);
  }
  initial_state.set_x(x);
  initial_state.set_y(y);
  initial_state.set_z(z);
  initial_state.set_vx(vx);
  initial_state.set_vy(vy);
  initial_state.set_vz(vz);
  CHECK_NOTNULL(plugin)->
      DirectlyInsertCelestial(
          celestial_index,
          parent_index,
          SolarSystem<Barycentric>::MakeDegreesOfFreedom(initial_state),
          SolarSystem<Barycentric>::MakeMassiveBody(gravity_model));
  return m.Return();
}

void principia__InsertCelestial(Plugin* const plugin,
                                int const celestial_index,
                                double const gravitational_parameter,
                                int const parent_index,
                                QP const from_parent) {
  journal::Method<journal::InsertCelestial> m({plugin,
                                               celestial_index,
                                               gravitational_parameter,
                                               parent_index,
                                               from_parent});
  CHECK_NOTNULL(plugin)->InsertCelestial(
      celestial_index,
      gravitational_parameter * SIUnit<GravitationalParameter>(),
      parent_index,
      RelativeDegreesOfFreedom<AliceSun>(
          Displacement<AliceSun>(ToR3Element(from_parent.q) * Metre),
          Velocity<AliceSun>(ToR3Element(from_parent.p) * (Metre / Second))));
  return m.Return();
}

void principia__InsertSun(Plugin* const plugin,
                          int const celestial_index,
                          double const gravitational_parameter) {
  journal::Method<journal::InsertSun> m({plugin,
                                         celestial_index,
                                         gravitational_parameter});
  CHECK_NOTNULL(plugin)->InsertSun(
      celestial_index,
      gravitational_parameter * SIUnit<GravitationalParameter>());
  return m.Return();
}

void principia__UpdateCelestialHierarchy(Plugin const* const plugin,
                                         int const celestial_index,
                                         int const parent_index) {
  journal::Method<journal::UpdateCelestialHierarchy> m({plugin,
                                                        celestial_index,
                                                        parent_index});
  CHECK_NOTNULL(plugin)->UpdateCelestialHierarchy(celestial_index,
                                                  parent_index);
  return m.Return();
}

void principia__EndInitialization(Plugin* const plugin) {
  journal::Method<journal::EndInitialization> m({plugin});
  CHECK_NOTNULL(plugin)->EndInitialization();
  return m.Return();
}

bool principia__InsertOrKeepVessel(Plugin* const plugin,
                                   char const* const vessel_guid,
                                   int const parent_index) {
  journal::Method<journal::InsertOrKeepVessel> m({plugin,
                                                  vessel_guid,
                                                  parent_index});
  return m.Return(
      CHECK_NOTNULL(plugin)->InsertOrKeepVessel(vessel_guid, parent_index));
}

void principia__SetVesselStateOffset(Plugin* const plugin,
                                     char const* const vessel_guid,
                                     QP const from_parent) {
  journal::Method<journal::SetVesselStateOffset> m({plugin,
                                                    vessel_guid,
                                                    from_parent});
  CHECK_NOTNULL(plugin)->SetVesselStateOffset(
      vessel_guid,
      RelativeDegreesOfFreedom<AliceSun>(
          Displacement<AliceSun>(ToR3Element(from_parent.q) * Metre),
          Velocity<AliceSun>(ToR3Element(from_parent.p) * (Metre / Second))));
  return m.Return();
}

void principia__AdvanceTime(Plugin* const plugin,
                            double const t,
                            double const planetarium_rotation) {
  journal::Method<journal::AdvanceTime> m({plugin, t, planetarium_rotation});
  Instant const t0;
  CHECK_NOTNULL(plugin)->AdvanceTime(t0 + t * Second,
                                     planetarium_rotation * Degree);
  return m.Return();
}

void principia__ForgetAllHistoriesBefore(Plugin* const plugin,
                                         double const t) {
  journal::Method<journal::ForgetAllHistoriesBefore> m({plugin, t});
  Instant const t0;
  CHECK_NOTNULL(plugin)->ForgetAllHistoriesBefore(t0 + t * Second);
  return m.Return();
}

QP principia__VesselFromParent(Plugin const* const plugin,
                               char const* const vessel_guid) {
  journal::Method<journal::VesselFromParent> m({plugin, vessel_guid});
  RelativeDegreesOfFreedom<AliceSun> const result =
      CHECK_NOTNULL(plugin)->VesselFromParent(vessel_guid);
  return m.Return({ToXYZ(result.displacement().coordinates() / Metre),
                   ToXYZ(result.velocity().coordinates() / (Metre / Second))});
}

QP principia__CelestialFromParent(Plugin const* const plugin,
                                  int const celestial_index) {
  journal::Method<journal::CelestialFromParent> m({plugin, celestial_index});
  RelativeDegreesOfFreedom<AliceSun> const result =
      CHECK_NOTNULL(plugin)->CelestialFromParent(celestial_index);
  return m.Return({ToXYZ(result.displacement().coordinates() / Metre),
                   ToXYZ(result.velocity().coordinates() / (Metre / Second))});
}

NavigationFrame* principia__NewBodyCentredNonRotatingNavigationFrame(
    Plugin const* const plugin,
    int const reference_body_index) {
  journal::Method<journal::NewBodyCentredNonRotatingNavigationFrame> m(
      {plugin, reference_body_index});
  return m.Return(CHECK_NOTNULL(plugin)->
      NewBodyCentredNonRotatingNavigationFrame(reference_body_index).release());
}

NavigationFrame* principia__NewBarycentricRotatingNavigationFrame(
    Plugin const* const plugin,
    int const primary_index,
    int const secondary_index) {
  journal::Method<journal::NewBarycentricRotatingNavigationFrame> m(
      {plugin, primary_index, secondary_index});
  return m.Return(CHECK_NOTNULL(plugin)->
      NewBarycentricRotatingNavigationFrame(primary_index,
                                           secondary_index).release());
}

void principia__SetPlottingFrame(Plugin* const plugin,
                                 NavigationFrame** const navigation_frame) {
  journal::Method<journal::SetPlottingFrame> m({plugin, navigation_frame},
                                               {navigation_frame});
  CHECK_NOTNULL(plugin)->SetPlottingFrame(TakeOwnership(navigation_frame));
  return m.Return();
}

void principia__UpdatePrediction(Plugin const* const plugin,
                                 char const* const vessel_guid) {
  journal::Method<journal::UpdatePrediction> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin)->UpdatePrediction(vessel_guid);
  return m.Return();
}

LineAndIterator* principia__RenderedVesselTrajectory(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position) {
  journal::Method<journal::RenderedVesselTrajectory> m({plugin,
                                                        vessel_guid,
                                                        sun_world_position});
  RenderedTrajectory<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedVesselTrajectory(
          vessel_guid,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return m.Return(result.release());
}

bool principia__HasPrediction(Plugin const* const plugin,
                              char const* const vessel_guid) {
  journal::Method<journal::HasPrediction> m({plugin, vessel_guid});
  return m.Return(CHECK_NOTNULL(plugin)->HasPrediction(vessel_guid));
}

LineAndIterator* principia__RenderedPrediction(
    Plugin* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position) {
  journal::Method<journal::RenderedPrediction> m({plugin,
                                                  vessel_guid,
                                                  sun_world_position});
  RenderedTrajectory<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedPrediction(
          vessel_guid,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return m.Return(result.release());
}

int principia__FlightPlanSize(Plugin const* const plugin,
                              char const* const vessel_guid) {
  journal::Method<journal::FlightPlanSize> m({plugin, vessel_guid});
  return m.Return(CHECK_NOTNULL(plugin)->FlightPlanSize(vessel_guid));
}

LineAndIterator* principia__RenderedFlightPlan(
    Plugin* const plugin,
    char const* const vessel_guid,
    int const plan_phase,
    XYZ const sun_world_position) {
  journal::Method<journal::RenderedFlightPlan> m({plugin,
                                                  vessel_guid,
                                                  plan_phase,
                                                  sun_world_position});
  RenderedTrajectory<World> rendered_trajectory =
      CHECK_NOTNULL(plugin)->RenderedFlightPlan(
          vessel_guid,
          plan_phase,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return m.Return(result.release());
}

void principia__SetPredictionLength(Plugin* const plugin,
                                    double const t) {
  journal::Method<journal::SetPredictionLength> m({plugin, t});
  CHECK_NOTNULL(plugin)->SetPredictionLength(t * Second);
  return m.Return();
}

void principia__SetPredictionLengthTolerance(Plugin* const plugin,
                                             double const l) {
  journal::Method<journal::SetPredictionLengthTolerance> m({plugin, l});
  CHECK_NOTNULL(plugin)->SetPredictionLengthTolerance(l * Metre);
  return m.Return();
}

void principia__SetPredictionSpeedTolerance(Plugin* const plugin,
                                            double const v) {
  journal::Method<journal::SetPredictionSpeedTolerance> m({plugin, v});
  CHECK_NOTNULL(plugin)->SetPredictionSpeedTolerance(v * Metre / Second);
  return m.Return();
}

bool principia__HasVessel(Plugin* const plugin,
                          char const* const vessel_guid) {
  journal::Method<journal::HasVessel> m({plugin,  vessel_guid});
  return m.Return(CHECK_NOTNULL(plugin)->HasVessel(vessel_guid));
}

int principia__NumberOfSegments(
    LineAndIterator const* const line_and_iterator) {
  journal::Method<journal::NumberOfSegments> m({line_and_iterator});
  return m.Return(CHECK_NOTNULL(line_and_iterator)->rendered_trajectory.size());
}

XYZSegment principia__FetchAndIncrement(
    LineAndIterator* const line_and_iterator) {
  journal::Method<journal::FetchAndIncrement> m({line_and_iterator});
  CHECK_NOTNULL(line_and_iterator);
  CHECK(line_and_iterator->it != line_and_iterator->rendered_trajectory.end());
  LineSegment<World> const result = *line_and_iterator->it;
  ++line_and_iterator->it;
  return m.Return({ToXYZ((result.begin - World::origin).coordinates() / Metre),
                   ToXYZ((result.end - World::origin).coordinates() / Metre)});
}

bool principia__AtEnd(LineAndIterator const* const line_and_iterator) {
  journal::Method<journal::AtEnd> m({line_and_iterator});
  CHECK_NOTNULL(line_and_iterator);
  return m.Return(line_and_iterator->it ==
                  line_and_iterator->rendered_trajectory.end());
}

void principia__DeleteLineAndIterator(
    LineAndIterator** const line_and_iterator) {
  journal::Method<journal::DeleteLineAndIterator> m({line_and_iterator},
                                                    {line_and_iterator});
  TakeOwnership(line_and_iterator);
  return m.Return();
}

void principia__AddVesselToNextPhysicsBubble(Plugin* const plugin,
                                             char const* const vessel_guid,
                                             KSPPart const* const parts,
                                             int count) {
  journal::Method<journal::AddVesselToNextPhysicsBubble> m({plugin,
                                                            vessel_guid,
                                                            parts,
                                                            count});
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
  return m.Return();
}

bool principia__PhysicsBubbleIsEmpty(Plugin const* const plugin) {
  journal::Method<journal::PhysicsBubbleIsEmpty> m({plugin});
  return m.Return(CHECK_NOTNULL(plugin)->PhysicsBubbleIsEmpty());
}

XYZ principia__BubbleDisplacementCorrection(Plugin const* const plugin,
                                            XYZ const sun_position) {
  journal::Method<journal::BubbleDisplacementCorrection> m({plugin,
                                                            sun_position});
  Displacement<World> const result =
      CHECK_NOTNULL(plugin)->BubbleDisplacementCorrection(
          World::origin + Displacement<World>(
                              ToR3Element(sun_position) * Metre));
  return m.Return(ToXYZ(result.coordinates() / Metre));
}

XYZ principia__BubbleVelocityCorrection(Plugin const* const plugin,
                                        int const reference_body_index) {
  journal::Method<journal::BubbleVelocityCorrection> m({plugin,
                                                        reference_body_index});
  Velocity<World> const result =
      CHECK_NOTNULL(plugin)->BubbleVelocityCorrection(reference_body_index);
  return m.Return(ToXYZ(result.coordinates() / (Metre / Second)));
}

WXYZ principia__NavballOrientation(
    Plugin const* const plugin,
    XYZ const sun_world_position,
    XYZ const ship_world_position) {
  journal::Method<journal::NavballOrientation> m({plugin,
                                                  sun_world_position,
                                                  ship_world_position});
  FrameField<World> const frame_field = CHECK_NOTNULL(plugin)->Navball(
      World::origin +
          Displacement<World>(ToR3Element(sun_world_position) * Metre));
  return m.Return(ToWXYZ(
      frame_field(
          World::origin +
              Displacement<World>(
                  ToR3Element(ship_world_position) * Metre)).quaternion()));
}

XYZ principia__VesselTangent(Plugin const* const plugin,
                             char const* const vessel_guid) {
  journal::Method<journal::VesselTangent> m({plugin, vessel_guid});
  return m.Return(
      ToXYZ(CHECK_NOTNULL(plugin)->VesselTangent(vessel_guid).coordinates()));
}

XYZ principia__VesselNormal(Plugin const* const plugin,
                            char const* const vessel_guid) {
  journal::Method<journal::VesselNormal> m({plugin, vessel_guid});
  return m.Return(
      ToXYZ(CHECK_NOTNULL(plugin)->VesselNormal(vessel_guid).coordinates()));
}

XYZ principia__VesselBinormal(Plugin const* const plugin,
                              char const* const vessel_guid) {
  journal::Method<journal::VesselBinormal> m({plugin, vessel_guid});
  return m.Return(
      ToXYZ(CHECK_NOTNULL(plugin)->VesselBinormal(vessel_guid).coordinates()));
}

double principia__CurrentTime(Plugin const* const plugin) {
  journal::Method<journal::CurrentTime> m({plugin});
  return m.Return((CHECK_NOTNULL(plugin)->CurrentTime() - Instant()) / Second);
}

char const* principia__SerializePlugin(Plugin const* const plugin,
                                       PullSerializer** const serializer) {
  journal::Method<journal::SerializePlugin> m({plugin, serializer},
                                              {serializer});
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
    return m.Return(nullptr);
  }

  // Convert to hexadecimal and return to the client.
  std::int64_t const hexadecimal_size = (bytes.size << 1) + 1;
  UniqueBytes hexadecimal(hexadecimal_size);
  HexadecimalEncode(bytes, hexadecimal.get());
  hexadecimal.data.get()[hexadecimal_size - 1] = '\0';
  return m.Return(reinterpret_cast<char const*>(hexadecimal.data.release()));
}

void principia__DeletePluginSerialization(char const** const serialization) {
  journal::Method<journal::DeletePluginSerialization> m({serialization},
                                                        {serialization});
  LOG(INFO) << __FUNCTION__;
  TakeOwnershipArray(reinterpret_cast<uint8_t const**>(serialization));
  return m.Return();
}

void principia__DeserializePlugin(char const* const serialization,
                                  int const serialization_size,
                                  PushDeserializer** const deserializer,
                                  Plugin const** const plugin) {
  journal::Method<journal::DeserializePlugin> m({serialization,
                                                 serialization_size,
                                                 deserializer,
                                                 plugin},
                                                {deserializer, plugin});
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
                        [bytes]() { delete[] bytes; });

  // If the data was empty, delete the deserializer.  This ensures that
  // |*plugin| is filled.
  if (byte_size == 0) {
    TakeOwnership(deserializer);
  }
  return m.Return();
}

char const* principia__SayHello() {
  journal::Method<journal::SayHello> m;
  return m.Return("Hello from native C++!");
}

}  // namespace ksp_plugin
}  // namespace principia
