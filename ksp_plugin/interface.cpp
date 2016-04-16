
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
#include "physics/kepler_orbit.hpp"
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
using geometry::RadiusLatitudeLongitude;
using ksp_plugin::AliceSun;
using ksp_plugin::Barycentric;
using ksp_plugin::Part;
using ksp_plugin::Positions;
using ksp_plugin::World;
using physics::KeplerianElements;
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

namespace interface {

namespace {

int const kChunkSize = 64 << 10;
int const kNumberOfChunks = 8;

base::not_null<std::unique_ptr<MassiveBody>> MakeMassiveBody(
    char const* const gravitational_parameter,
    char const* const mean_radius,
    char const* const axis_right_ascension,
    char const* const axis_declination,
    char const* const j2,
    char const* const reference_radius) {
  serialization::GravityModel::Body gravity_model;
  gravity_model.set_gravitational_parameter(gravitational_parameter);
  gravity_model.set_mean_radius(mean_radius);
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
  return SolarSystem<Barycentric>::MakeMassiveBody(gravity_model);
}

}  // namespace

// Sets stderr to log INFO, and redirects stderr, which Unity does not log, to
// "<KSP directory>/stderr.log".  This provides an easily accessible file
// containing a sufficiently verbose log of the latest session, instead of
// requiring users to dig in the archive of all past logs at all severities.
// This archive is written to
// "<KSP directory>/glog/Principia/<SEVERITY>.<date>-<time>.<pid>",
// where date and time are in ISO 8601 basic format.
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

// If |activate| is true and there is no active journal, create one and
// activate it.  If |activate| is false and there is an active journal,
// deactivate it.  Does nothing if there is already a journal in the desired
// state.  |verbose| causes methods to be output in the INFO log before being
// executed.
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

// Log messages at a level |<= max_severity| are buffered.
// Log messages at a higher level are flushed immediately.
void principia__SetBufferedLogging(int const max_severity) {
  journal::Method<journal::SetBufferedLogging> m({max_severity});
  FLAGS_logbuflevel = max_severity;
  return m.Return();
}

int principia__GetBufferedLogging() {
  journal::Method<journal::GetBufferedLogging> m;
  return m.Return(FLAGS_logbuflevel);
}

// Sets the maximum number of seconds which logs may be buffered for.
void principia__SetBufferDuration(int const seconds) {
  journal::Method<journal::SetBufferDuration> m({seconds});
  FLAGS_logbufsecs = seconds;
  return m.Return();
}

int principia__GetBufferDuration() {
  journal::Method<journal::GetBufferDuration> m;
  return m.Return(FLAGS_logbufsecs);
}

// Log suppression level: messages logged at a lower level than this are
// suppressed.
void principia__SetSuppressedLogging(int const min_severity) {
  journal::Method<journal::SetSuppressedLogging> m({min_severity});
  FLAGS_minloglevel = min_severity;
  return m.Return();
}

int principia__GetSuppressedLogging() {
  journal::Method<journal::GetSuppressedLogging> m;
  return m.Return(FLAGS_minloglevel);
}

// Show all VLOG(m) messages for |m <= level|.
void principia__SetVerboseLogging(int const level) {
  journal::Method<journal::SetVerboseLogging> m({level});
  FLAGS_v = level;
  return m.Return();
}

int principia__GetVerboseLogging() {
  journal::Method<journal::GetVerboseLogging> m;
  return m.Return(FLAGS_v);
}

// Make it so that all log messages of at least |min_severity| are logged to
// stderr (in addition to logging to the usual log file(s)).
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

// Exports |LOG(SEVERITY) << text| for fast logging from the C# adapter.
// This will always evaluate its argument even if the corresponding log severity
// is disabled, so it is less efficient than LOG(INFO).  It will not report the
// line and file of the caller.
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

// Returns a pointer to a plugin constructed with the arguments given.
// The caller takes ownership of the result.
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

// Deletes and nulls |*plugin|.
// |plugin| must not be null.  No transfer of ownership of |*plugin|, takes
// ownership of |**plugin|.
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

void principia__InsertCelestialAbsoluteCartesian(
    Plugin* const plugin,
    int const celestial_index,
    int const* const parent_index,
    char const* const gravitational_parameter,
    char const* const mean_radius,
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
  journal::Method<journal::InsertCelestialAbsoluteCartesian> m(
      {plugin,
       celestial_index,
       parent_index,
       gravitational_parameter,
       mean_radius,
       axis_right_ascension,
       axis_declination,
       j2,
       reference_radius,
       x, y, z,
       vx, vy, vz});
  serialization::InitialState::Body initial_state;
  initial_state.set_x(x);
  initial_state.set_y(y);
  initial_state.set_z(z);
  initial_state.set_vx(vx);
  initial_state.set_vy(vy);
  initial_state.set_vz(vz);
  CHECK_NOTNULL(plugin)
      ->InsertCelestialAbsoluteCartesian(
          celestial_index,
          parent_index == nullptr
              ? std::experimental::nullopt
              : std::experimental::make_optional(*parent_index),
          SolarSystem<Barycentric>::MakeDegreesOfFreedom(initial_state),
          MakeMassiveBody(gravitational_parameter,
                          mean_radius,
                          axis_right_ascension,
                          axis_declination,
                          j2,
                          reference_radius));
  return m.Return();
}

void principia__InsertCelestialJacobiKeplerian(
    Plugin* const plugin,
    int const celestial_index,
    int const parent_index,
    char const* const gravitational_parameter,
    char const* const mean_radius,
    char const* const axis_right_ascension,
    char const* const axis_declination,
    char const* const j2,
    char const* const reference_radius,
    double const eccentricity,
    char const* const mean_motion,
    char const* const inclination,
    char const* const longitude_of_ascending_node,
    char const* const argument_of_periapsis,
    char const* const mean_anomaly) {
  journal::Method<journal::InsertCelestialJacobiKeplerian> m(
      {plugin,
       celestial_index,
       parent_index,
       gravitational_parameter,
       mean_radius,
       axis_right_ascension,
       axis_declination,
       j2,
       reference_radius,
       eccentricity,
       mean_motion,
       inclination,
       longitude_of_ascending_node,
       argument_of_periapsis,
       mean_anomaly});
  KeplerianElements<Barycentric> keplerian_elements;
  keplerian_elements.eccentricity = eccentricity;
  keplerian_elements.mean_motion = ParseQuantity<AngularFrequency>(mean_motion);
  keplerian_elements.inclination = ParseQuantity<Angle>(inclination);
  keplerian_elements.longitude_of_ascending_node =
      ParseQuantity<Angle>(longitude_of_ascending_node);
  keplerian_elements.argument_of_periapsis =
      ParseQuantity<Angle>(argument_of_periapsis);
  keplerian_elements.mean_anomaly = ParseQuantity<Angle>(mean_anomaly);
  CHECK_NOTNULL(plugin)
      ->InsertCelestialJacobiKeplerian(
          celestial_index,
          parent_index,
          keplerian_elements,
          MakeMassiveBody(gravitational_parameter,
                          mean_radius,
                          axis_right_ascension,
                          axis_declination,
                          j2,
                          reference_radius));
  return m.Return();
}

void principia__InsertSun(Plugin* const plugin,
                          int const celestial_index,
                          double const gravitational_parameter,
                          double const mean_radius) {
  journal::Method<journal::InsertSun> m({plugin,
                                         celestial_index,
                                         gravitational_parameter,
                                         mean_radius});
  CHECK_NOTNULL(plugin)->InsertSun(
      celestial_index,
      gravitational_parameter * SIUnit<GravitationalParameter>(),
      mean_radius * Metre);
  return m.Return();
}

// Calls |plugin->UpdateCelestialHierarchy| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
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

// Calls |plugin->EndInitialization|.
// |plugin| must not be null.  No transfer of ownership.
void principia__EndInitialization(Plugin* const plugin) {
  journal::Method<journal::EndInitialization> m({plugin});
  CHECK_NOTNULL(plugin)->EndInitialization();
  return m.Return();
}

// Calls |plugin->InsertOrKeepVessel| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
bool principia__InsertOrKeepVessel(Plugin* const plugin,
                                   char const* const vessel_guid,
                                   int const parent_index) {
  journal::Method<journal::InsertOrKeepVessel> m({plugin,
                                                  vessel_guid,
                                                  parent_index});
  return m.Return(
      CHECK_NOTNULL(plugin)->InsertOrKeepVessel(vessel_guid, parent_index));
}

// Calls |plugin->SetVesselStateOffset| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
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

// Calls |plugin->VesselFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
QP principia__VesselFromParent(Plugin const* const plugin,
                               char const* const vessel_guid) {
  journal::Method<journal::VesselFromParent> m({plugin, vessel_guid});
  RelativeDegreesOfFreedom<AliceSun> const result =
      CHECK_NOTNULL(plugin)->VesselFromParent(vessel_guid);
  return m.Return({ToXYZ(result.displacement().coordinates() / Metre),
                   ToXYZ(result.velocity().coordinates() / (Metre / Second))});
}

// Calls |plugin->CelestialFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
QP principia__CelestialFromParent(Plugin const* const plugin,
                                  int const celestial_index) {
  journal::Method<journal::CelestialFromParent> m({plugin, celestial_index});
  RelativeDegreesOfFreedom<AliceSun> const result =
      CHECK_NOTNULL(plugin)->CelestialFromParent(celestial_index);
  return m.Return({ToXYZ(result.displacement().coordinates() / Metre),
                   ToXYZ(result.velocity().coordinates() / (Metre / Second))});
}

// Calls |plugin->NewBodyCentredNonRotatingFrame| with the arguments given.
// |plugin| must not be null.  The caller gets ownership of the returned object.
// TODO(phl): The parameter should be named |centre_index|.
NavigationFrame* principia__NewBodyCentredNonRotatingNavigationFrame(
    Plugin const* const plugin,
    int const reference_body_index) {
  journal::Method<journal::NewBodyCentredNonRotatingNavigationFrame> m(
      {plugin, reference_body_index});
  return m.Return(CHECK_NOTNULL(plugin)->
      NewBodyCentredNonRotatingNavigationFrame(reference_body_index).release());
}

// Calls |plugin->NewBarycentricRotatingFrame| with the arguments given.
// |plugin| must not be null.  The caller gets ownership of the returned object.
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

// Calls |plugin| to create a |NavigationFrame| using the given |parameters|.
NavigationFrame* principia__NewNavigationFrame(
    Plugin const* const plugin,
    NavigationFrameParameters const parameters) {
  journal::Method<journal::NewNavigationFrame> m({plugin, parameters});
  return m.Return(NewNavigationFrame(plugin, parameters).release());
}

// |navigation_frame| must not be null.  No transfer of ownership of
// |*navigation_frame|, takes ownership of |**navigation_frame|, nulls
// |*navigation_frame|.
void principia__SetPlottingFrame(Plugin* const plugin,
                                 NavigationFrame** const navigation_frame) {
  journal::Method<journal::SetPlottingFrame> m({plugin, navigation_frame},
                                               {navigation_frame});
  CHECK_NOTNULL(plugin)->SetPlottingFrame(TakeOwnership(navigation_frame));
  return m.Return();
}

// Returns the frame last set by |plugin->SetPlottingFrame|.  No transfer of
// ownership.  The returned pointer is never null.
NavigationFrame const* principia__GetPlottingFrame(Plugin const* const plugin) {
  journal::Method<journal::GetPlottingFrame> m({plugin});
  return m.Return(CHECK_NOTNULL(plugin)->GetPlottingFrame());
}

void principia__UpdatePrediction(Plugin const* const plugin,
                                 char const* const vessel_guid) {
  journal::Method<journal::UpdatePrediction> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin)->UpdatePrediction(vessel_guid);
  return m.Return();
}

// Returns the result of |plugin->RenderedVesselTrajectory| called with the
// arguments given, together with an iterator to its beginning.
// |plugin| must not be null.  No transfer of ownership of |plugin|.  The caller
// gets ownership of the result.  |frame| must not be null.  No transfer of
// ownership of |frame|.
Iterator* principia__RenderedVesselTrajectory(Plugin const* const plugin,
                                              char const* const vessel_guid,
                                              XYZ const sun_world_position) {
  journal::Method<journal::RenderedVesselTrajectory> m({plugin,
                                                        vessel_guid,
                                                        sun_world_position});
  Positions<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedVesselTrajectory(
          vessel_guid,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  return m.Return(new TypedIterator<Positions<World>>(
      std::move(rendered_trajectory)));
}

Iterator* principia__RenderedPrediction(Plugin* const plugin,
                                        char const* const vessel_guid,
                                        XYZ const sun_world_position) {
  journal::Method<journal::RenderedPrediction> m({plugin,
                                                  vessel_guid,
                                                  sun_world_position});
  Positions<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedPrediction(
          vessel_guid,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  return m.Return(new TypedIterator<Positions<World>>(
      std::move(rendered_trajectory)));
}

void principia__RenderedPredictionApsides(Plugin const* const plugin,
                                          char const* const vessel_guid,
                                          int const celestial_index,
                                          XYZ const sun_world_position,
                                          Iterator** const apoapsides,
                                          Iterator** const periapsides) {
  journal::Method<journal::RenderedPredictionApsides> m(
      {plugin, vessel_guid, celestial_index, sun_world_position},
      {apoapsides, periapsides});
  CHECK_NOTNULL(plugin);
  auto const& prediction = plugin->GetVessel(vessel_guid)->prediction();
  Position<World> q_sun =
      World::origin +
      Displacement<World>(ToR3Element(sun_world_position) * Metre);
  Positions<World> rendered_apoapsides;
  Positions<World> rendered_periapsides;
  plugin->ComputeAndRenderApsides(celestial_index,
                                  prediction.Fork(),
                                  prediction.End(),
                                  q_sun,
                                  rendered_apoapsides,
                                  rendered_periapsides);
  *apoapsides = new TypedIterator<Positions<World>>(
      std::move(rendered_apoapsides));
  *periapsides = new TypedIterator<Positions<World>>(
      std::move(rendered_periapsides));
  return m.Return();
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
                part->mass_in_tonnes * Tonne,
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

// |plugin| must not be null.  The caller takes ownership of the result, except
// when it is null (at the end of the stream).  No transfer of ownership of
// |*plugin|.  |*serializer| must be null on the first call and must be passed
// unchanged to the successive calls; its ownership is not transferred.
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

// Deletes and nulls |*serialization|.
// |serialization| must not be null.  No transfer of ownership of
// |*serialization|, takes ownership of |**serialization|.
void principia__DeletePluginSerialization(char const** const serialization) {
  journal::Method<journal::DeletePluginSerialization> m({serialization},
                                                        {serialization});
  LOG(INFO) << __FUNCTION__;
  TakeOwnershipArray(reinterpret_cast<uint8_t const**>(serialization));
  return m.Return();
}

// The caller takes ownership of |**plugin| when it is not null.  No transfer of
// ownership of |*serialization| or |**deserializer|.  |*deserializer| and
// |*plugin| must be null on the first call and must be passed unchanged to the
// successive calls.  The caller must perform an extra call with
// |serialization_size| set to 0 to indicate the end of the input stream.  When
// this last call returns, |*plugin| is not null and may be used by the caller.
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

// Says hello, convenient for checking that calls to the DLL work.
char const* principia__SayHello() {
  journal::Method<journal::SayHello> m;
  return m.Return("Hello from native C++!");
}

void principia__GetVersion(
    char const** const build_date,
    char const** const version) {
  journal::Method<journal::GetVersion> m({build_date, version});
  *CHECK_NOTNULL(build_date) = base::kBuildDate;
  *CHECK_NOTNULL(version) = base::kVersion;
  return m.Return();
}

}  // namespace interface
}  // namespace principia
