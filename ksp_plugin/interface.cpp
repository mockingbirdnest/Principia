#include "ksp_plugin/interface.hpp"

#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#define MICROSOFT_WINDOWS_WINBASE_H_DEFINE_INTERLOCKED_CPLUSPLUS_OVERLOADS 0
#if OS_WIN
#include <windows.h>
#include <psapi.h>
#endif

#include "absl/strings/str_split.h"
#include "astronomy/epoch.hpp"
#include "astronomy/time_scales.hpp"
#include "base/array.hpp"
#include "base/base64.hpp"
#include "base/cpuid.hpp"
#include "base/encoder.hpp"
#include "base/fingerprint2011.hpp"
#include "base/flags.hpp"
#include "base/hexadecimal.hpp"
#include "base/macros.hpp"  // 🧙 For NAMED.
#include "base/not_null.hpp"
#include "base/optional_logging.hpp"  // 🧙 For logging.
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "base/serialization.hpp"
#include "base/version.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "gipfeli/gipfeli.h"
#include "google/protobuf/arena.h"
#include "integrators/integrators.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "journal/recorder.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/part.hpp"
#include "numerics/elementary_functions.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/ephemeris.hpp"
#include "physics/frame_field.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rotating_body.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "physics/tensors.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/astronomy.pb.h"
#include "serialization/geometry.pb.h"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace interface {

using ::google::protobuf::Arena;
using ::google::protobuf::ArenaOptions;
using ::operator<<;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_array;
using namespace principia::base::_base64;
using namespace principia::base::_cpuid;
using namespace principia::base::_encoder;
using namespace principia::base::_fingerprint2011;
using namespace principia::base::_flags;
using namespace principia::base::_hexadecimal;
using namespace principia::base::_not_null;
using namespace principia::base::_pull_serializer;
using namespace principia::base::_push_deserializer;
using namespace principia::base::_serialization;
using namespace principia::base::_version;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_rotation;
using namespace principia::integrators::_integrators;
using namespace principia::journal::_method;
using namespace principia::journal::_recorder;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_identification;
using namespace principia::ksp_plugin::_iterators;
using namespace principia::ksp_plugin::_part;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_segment;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_frame_field;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_oblate_body;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::physics::_tensors;
using namespace principia::quantities::_astronomy;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_parser;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

namespace {

constexpr char gipfeli_compressor[] = "gipfeli";

constexpr char base64_encoder[] = "base64";
constexpr char hexadecimal_encoder[] = "hexadecimal";

constexpr int chunk_size = 64 << 10;
constexpr int number_of_chunks = 8;

not_null<Arena*> arena = []() {
  ArenaOptions options;
  options.initial_block_size = chunk_size;
  options.max_block_size = 16 * chunk_size;
  return new Arena(options);
}();

Ephemeris<Barycentric>::AccuracyParameters MakeAccuracyParameters(
    ConfigurationAccuracyParameters const& parameters) {
  return Ephemeris<Barycentric>::AccuracyParameters(
      ParseQuantity<Length>(parameters.fitting_tolerance),
      ParseQuantity<double>(parameters.geopotential_tolerance));
}

Ephemeris<Barycentric>::AdaptiveStepParameters MakeAdaptiveStepParameters(
    ConfigurationAdaptiveStepParameters const& parameters) {
  // It is erroneous for a psychohistory integration to fail, so the `max_steps`
  // must be unlimited.
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      ParseAdaptiveStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation>(
          parameters.adaptive_step_size_integrator),
      /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
      ParseQuantity<Length>(parameters.length_integration_tolerance),
      ParseQuantity<Speed>(parameters.speed_integration_tolerance));
}

DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
MakeDownsamplingParameters(
    ConfigurationDownsamplingParameters const& parameters) {
  return DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters{
      std::stoi(parameters.max_dense_intervals),
      ParseQuantity<Length>(parameters.tolerance)};
}

Ephemeris<Barycentric>::FixedStepParameters MakeFixedStepParameters(
    ConfigurationFixedStepParameters const& parameters) {
  return Ephemeris<Barycentric>::FixedStepParameters(
      ParseFixedStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation>(
          parameters.fixed_step_size_integrator),
      ParseQuantity<Time>(parameters.integration_step_size));
}

serialization::OblateBody::Geopotential MakeGeopotential(
    BodyGeopotentialElement const* const* const geopotential) {
  // Make sure that we generate at most one row per degree.
  std::map<int, serialization::OblateBody::Geopotential::GeopotentialRow> rows;
  for (BodyGeopotentialElement const* const* e = geopotential;
       *e != nullptr;
       ++e) {
    BodyGeopotentialElement const& element = **e;
    serialization::OblateBody::Geopotential::GeopotentialRow::GeopotentialColumn
        column;
    int const degree = std::stoi(element.degree);
    int const order = std::stoi(element.order);
    column.set_order(order);
    if (element.j != nullptr) {
      CHECK(element.cos == nullptr);
      double const j = ParseQuantity<double>(element.j);
      column.set_j(j);
    }
    if (element.cos != nullptr) {
      CHECK(element.j == nullptr);
      double const cos = ParseQuantity<double>(element.cos);
      column.set_cos(cos);
    }
    double const sin = ParseQuantity<double>(element.sin);
    column.set_sin(sin);
    *rows[degree].add_column() = column;
    rows[degree].set_degree(degree);
  }

  serialization::OblateBody::Geopotential result;
  for (auto const& [_, row] : rows) {
    *result.add_row() = row;
  }
  return result;
}

serialization::GravityModel::Body MakeGravityModel(
    BodyParameters const& body_parameters) {
  // Logging operators would dereference a null C string.
  auto const make_optional_c_string = [](char const* const c_string) {
    return (c_string == nullptr) ? std::nullopt
                                 : std::make_optional(c_string);
  };
  auto const make_optional_geopotential_size =
      [](BodyGeopotentialElement const* const* const geopotential)
      -> std::optional<std::int64_t> {
    if (geopotential == nullptr) {
      return std::nullopt;
    } else {
      std::int64_t size = 0;
      for (BodyGeopotentialElement const* const* e = geopotential;
           *e != nullptr;
           ++e) {
        ++size;
      }
      return size;
    }
  };
  LOG(INFO)
      << __FUNCTION__ << "\n"
      << NAMED(make_optional_c_string(body_parameters.gravitational_parameter))
      << "\n"
      << NAMED(body_parameters.reference_instant) << "\n"
      << NAMED(make_optional_c_string(body_parameters.min_radius)) << "\n"
      << NAMED(make_optional_c_string(body_parameters.mean_radius)) << "\n"
      << NAMED(make_optional_c_string(body_parameters.max_radius)) << "\n"
      << NAMED(make_optional_c_string(body_parameters.axis_right_ascension))
      << "\n"
      << NAMED(make_optional_c_string(body_parameters.axis_declination)) << "\n"
      << NAMED(make_optional_c_string(body_parameters.reference_angle)) << "\n"
      << NAMED(make_optional_c_string(body_parameters.angular_frequency))
      << "\n"
      << NAMED(make_optional_c_string(body_parameters.j2)) << "\n"
      << NAMED(make_optional_c_string(body_parameters.reference_radius)) << "\n"
      << NAMED(make_optional_geopotential_size(body_parameters.geopotential));
  serialization::GravityModel::Body gravity_model;
  gravity_model.set_name(body_parameters.name);
  gravity_model.set_gravitational_parameter(
      body_parameters.gravitational_parameter);
  if (body_parameters.reference_instant != nullptr) {
    gravity_model.set_reference_instant(body_parameters.reference_instant);
  }
  if (body_parameters.min_radius != nullptr) {
    gravity_model.set_min_radius(body_parameters.min_radius);
  }
  if (body_parameters.mean_radius != nullptr) {
    gravity_model.set_mean_radius(body_parameters.mean_radius);
  }
  if (body_parameters.max_radius != nullptr) {
    gravity_model.set_max_radius(body_parameters.max_radius);
  }
  if (body_parameters.axis_right_ascension != nullptr) {
    gravity_model.set_axis_right_ascension(
        body_parameters.axis_right_ascension);
  }
  if (body_parameters.axis_declination != nullptr) {
    gravity_model.set_axis_declination(body_parameters.axis_declination);
  }
  if (body_parameters.reference_angle != nullptr) {
    gravity_model.set_reference_angle(body_parameters.reference_angle);
  }
  if (body_parameters.angular_frequency!= nullptr) {
    gravity_model.set_angular_frequency(
        body_parameters.angular_frequency);
  }
  if (body_parameters.reference_radius != nullptr) {
    gravity_model.set_reference_radius(body_parameters.reference_radius);
  }
  if (body_parameters.j2 != nullptr) {
    gravity_model.set_j2(ParseQuantity<double>(body_parameters.j2));
  }
  if (body_parameters.geopotential != nullptr) {
    *gravity_model.mutable_geopotential() =
        MakeGeopotential(body_parameters.geopotential);
  }
  LOG(INFO) << "Fingerprint " << std::setw(16) << std::hex << std::uppercase
            << Fingerprint2011(SerializeAsBytes(gravity_model).get())
            << " for " << gravity_model.name();
  return gravity_model;
}

std::unique_ptr<google::compression::Compressor> NewCompressor(
    std::string_view const compressor) {
  if (compressor.empty()) {
    return nullptr;
  } else if (compressor == gipfeli_compressor) {
    return google::compression::NewGipfeliCompressor();
  } else {
    LOG(FATAL) << "Unknown compressor " << compressor;
  }
}

Encoder<char, /*null_terminated=*/true>*
NewEncoder(std::string_view const encoder) {
  if (encoder == hexadecimal_encoder) {
    static auto* const encoder =
        new HexadecimalEncoder</*null_terminated=*/true>;
    return encoder;
  } else if (encoder == base64_encoder) {
    static auto* const encoder =
        new Base64Encoder</*null_terminated=*/true>;
    return encoder;
  } else {
    LOG(FATAL) << "Unknown encoder " << encoder;
  }
}

}  // namespace

void __cdecl principia__ActivatePlayer() {
  Vessel::MakeSynchronous();
}

// If `activate` is true and there is no active journal, create one and
// activate it.  If `activate` is false and there is an active journal,
// deactivate it.  Does nothing if there is already a journal in the desired
// state.  `verbose` causes methods to be output in the INFO log before being
// executed.
void __cdecl principia__ActivateRecorder(bool const activate) {
  // NOTE: Do not journal!  You'd end up with half a message in the journal and
  // that would cause trouble.
  if (activate && !Recorder::IsActivated()) {
    // Build a name somewhat similar to that of the log files.
    auto const now = std::chrono::system_clock::now();
    std::time_t const time = std::chrono::system_clock::to_time_t(now);
    std::tm* const localtime = std::localtime(&time);
    std::stringstream name;
    name << std::put_time(localtime, "JOURNAL.%Y%m%d-%H%M%S");
    Recorder* const recorder = new Recorder(
        std::filesystem::path("glog") / "Principia" / name.str());
    Vessel::MakeSynchronous();
    Recorder::Activate(recorder);
  } else if (!activate && Recorder::IsActivated()) {
    Recorder::Deactivate();
    Vessel::MakeAsynchronous();
  }
}

XYZ __cdecl principia__AngularMomentumFromAngularVelocity(
    XYZ world_angular_velocity,
    XYZ moments_of_inertia_in_tonnes,
    WXYZ principal_axes_rotation,
    WXYZ part_rotation) {
  journal::Method<journal::AngularMomentumFromAngularVelocity> m(
      {world_angular_velocity,
       moments_of_inertia_in_tonnes,
       principal_axes_rotation,
       part_rotation});
  using PartPrincipalAxes = Frame<serialization::Frame::PhysicsTag,
                                  Arbitrary,
                                  Handedness::Left,
                                  serialization::Frame::PRINCIPAL_AXES>;

  auto const angular_velocity =
      FromXYZ<AngularVelocity<World>>(world_angular_velocity);

  static constexpr MomentOfInertia zero;
  auto const moments_of_inertia =
      FromXYZ<R3Element<MomentOfInertia>>({moments_of_inertia_in_tonnes.x,
                                           moments_of_inertia_in_tonnes.y,
                                           moments_of_inertia_in_tonnes.z});
  InertiaTensor<PartPrincipalAxes> const inertia_tensor_in_princial_axes(
      R3x3Matrix<MomentOfInertia>({moments_of_inertia.x, zero, zero},
                                  {zero, moments_of_inertia.y, zero},
                                  {zero, zero, moments_of_inertia.z}));

  Rotation<PartPrincipalAxes, RigidPart> const principal_axes_to_part(
      FromWXYZ(principal_axes_rotation));
  Rotation<RigidPart, World> const part_to_world(FromWXYZ(part_rotation));

  InertiaTensor<World> const inertia_tensor =
      part_to_world(principal_axes_to_part(inertia_tensor_in_princial_axes));

  Bivector<AngularMomentum, World> angular_momentum =
      inertia_tensor * angular_velocity;

  return m.Return(ToXYZ(angular_momentum));
}

void __cdecl principia__AdvanceTime(Plugin* const plugin,
                                    double const t,
                                    double const planetarium_rotation) {
  journal::Method<journal::AdvanceTime> m({plugin, t, planetarium_rotation});
  CHECK_NOTNULL(plugin);
  plugin->AdvanceTime(FromGameTime(*plugin, t), planetarium_rotation * Degree);
  return m.Return();
}

void __cdecl principia__CatchUpLaggingVessels(
    Plugin* const plugin,
    Iterator** const collided_vessels) {
  journal::Method<journal::CatchUpLaggingVessels> m({plugin},
                                                    {collided_vessels});
  CHECK_NOTNULL(plugin);
  VesselSet collided_vessel_set;
  plugin->CatchUpLaggingVessels(collided_vessel_set);
  *collided_vessels =
      new TypedIterator<VesselSet>(std::move(collided_vessel_set), plugin);
  return m.Return();
}

// Calls `plugin->CelestialFromParent` with the arguments given.
// `plugin` must not be null.  No transfer of ownership.
QP __cdecl principia__CelestialFromParent(Plugin const* const plugin,
                                          int const celestial_index) {
  journal::Method<journal::CelestialFromParent> m({plugin, celestial_index});
  CHECK_NOTNULL(plugin);
  return m.Return(ToQP(plugin->CelestialFromParent(celestial_index)));
}

double __cdecl principia__CelestialInitialRotationInDegrees(
    Plugin const* const plugin,
    int const celestial_index) {
  journal::Method<journal::CelestialInitialRotationInDegrees> m(
      {plugin, celestial_index});
  CHECK_NOTNULL(plugin);
  return m.Return(plugin->CelestialInitialRotation(celestial_index) / Degree);
}

WXYZ __cdecl principia__CelestialRotation(Plugin const* const plugin,
                                          int const index) {
  journal::Method<journal::CelestialRotation> m({plugin, index});
  CHECK_NOTNULL(plugin);
  return m.Return(ToWXYZ(plugin->CelestialRotation(index).quaternion()));
}

double __cdecl principia__CelestialRotationPeriod(
    Plugin const* const plugin,
    int const celestial_index) {
  journal::Method<journal::CelestialRotationPeriod> m(
      {plugin, celestial_index});
  CHECK_NOTNULL(plugin);
  return m.Return(plugin->CelestialRotationPeriod(celestial_index) / Second);
}

WXYZ __cdecl principia__CelestialSphereRotation(Plugin const* const plugin) {
  journal::Method<journal::CelestialSphereRotation> m({plugin});
  CHECK_NOTNULL(plugin);
  return m.Return(ToWXYZ(plugin->CelestialSphereRotation().quaternion()));
}

QP __cdecl principia__CelestialWorldDegreesOfFreedom(Plugin const* const plugin,
                                                     int const index,
                                                     Origin const origin,
                                                     double const time) {
  journal::Method<journal::CelestialWorldDegreesOfFreedom> m(
      {plugin, index, origin, time});
  CHECK_NOTNULL(plugin);
  return m.Return(ToQP(
      plugin->CelestialWorldDegreesOfFreedom(
          index,
          plugin->BarycentricToWorld(
              origin.reference_part_is_unmoving,
              origin.reference_part_id,
              origin.reference_part_is_at_origin
                  ? std::nullopt
                  : std::make_optional(
                        FromXYZ<Position<World>>(
                            origin.main_body_centre_in_world))),
          FromGameTime(*plugin, time))));
}

void __cdecl principia__ClearFlags() {
  journal::Method<journal::ClearFlags> m;
  Flags::Clear();
  return m.Return();
}

void __cdecl principia__ClearWorldRotationalReferenceFrame(
    Plugin* const plugin) {
  journal::Method<journal::ClearWorldRotationalReferenceFrame> m({plugin});
  CHECK_NOTNULL(plugin);
  plugin->ClearWorldRotationalReferenceFrame();
  return m.Return();
}

double __cdecl principia__CurrentTime(Plugin const* const plugin) {
  journal::Method<journal::CurrentTime> m({plugin});
  CHECK_NOTNULL(plugin);
  return m.Return(ToGameTime(*plugin, plugin->CurrentTime()));
}

void __cdecl principia__DeleteInterchange(void const** const native_pointer) {
  journal::Method<journal::DeleteInterchange> m({native_pointer},
                                                {native_pointer});
  CHECK_NOTNULL(native_pointer);
  ::operator delete(const_cast<void*>(*native_pointer));
  *native_pointer = nullptr;
  return m.Return();
}

// Deletes and nulls `*plugin`.
// `plugin` must not be null.  No transfer of ownership of `*plugin`, takes
// ownership of `**plugin`.
void __cdecl principia__DeletePlugin(Plugin const** const plugin) {
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

// Deletes and nulls `*native_string`.  `native_string` must not be null.  No
// transfer of ownership of `*native_string`, takes ownership of
// `**native_string`.
void __cdecl principia__DeleteString(char const** const native_string) {
  journal::Method<journal::DeleteString> m({native_string}, {native_string});
  TakeOwnershipArray(native_string);
  return m.Return();
}

// Same as above, but for char16_t.
void __cdecl principia__DeleteU16String(char16_t const** const native_string) {
  journal::Method<journal::DeleteU16String> m({native_string},
                                              {native_string});
  TakeOwnershipArray(native_string);
  return m.Return();
}

// The caller takes ownership of `**plugin` when it is not null.  No transfer of
// ownership of `*serialization` or `**deserializer`.  `*deserializer` and
// `*plugin` must be null on the first call and must be passed unchanged to the
// successive calls.  The caller must perform an extra call with
// `serialization_size` set to 0 to indicate the end of the input stream.  When
// this last call returns, `*plugin` is not null and may be used by the caller.
void __cdecl principia__DeserializePlugin(
    char const* const serialization,
    PushDeserializer** const deserializer,
    Plugin const** const plugin,
    char const* const compressor,
    char const* const encoder) {
  journal::Method<journal::DeserializePlugin> m({serialization,
                                                 deserializer,
                                                 plugin,
                                                 compressor,
                                                 encoder},
                                                {deserializer,
                                                 plugin});
  CHECK_NOTNULL(serialization);
  CHECK_NOTNULL(deserializer);
  CHECK_NOTNULL(plugin);

  // Create and start a deserializer if the caller didn't provide one.
  if (*deserializer == nullptr) {
    LOG(INFO) << "Begin plugin deserialization";
    *deserializer = new PushDeserializer(chunk_size,
                                         number_of_chunks,
                                         NewCompressor(compressor));
    CHECK_NOTNULL(arena);
    not_null<serialization::Plugin*> const message =
        Arena::CreateMessage<serialization::Plugin>(arena);
    (*deserializer)->Start(
        message,
        [plugin](google::protobuf::Message const& message) {
          *plugin = Plugin::ReadFromMessage(
              static_cast<serialization::Plugin const&>(message)).release();
        });
  }

  // Decode the representation.
  auto bytes = NewEncoder(encoder)->Decode({serialization,
                                            std::strlen(serialization)});
  auto const bytes_size = bytes.size;
  (*deserializer)->Push(std::move(bytes));

  // If the data was empty, delete the deserializer.  This ensures that
  // `*plugin` is filled.
  if (bytes_size == 0) {
    LOG(INFO) << "End plugin deserialization";
    TakeOwnership(deserializer);
    arena->Reset();
  }
  return m.Return();
}

// Calls `plugin->EndInitialization`.
// `plugin` must not be null.  No transfer of ownership.
void __cdecl principia__EndInitialization(Plugin* const plugin) {
  journal::Method<journal::EndInitialization> m({plugin});
  CHECK_NOTNULL(plugin);
  plugin->EndInitialization();
  return m.Return();
}

int __cdecl principia__EquipotentialCount(Plugin* const plugin) {
  journal::Method<journal::EquipotentialCount> m({plugin});
  CHECK_NOTNULL(plugin);
  plugin->geometric_potential_plotter().RefreshEquipotentials();
  auto const* equipotentials =
      plugin->geometric_potential_plotter().equipotentials();
  PlottingFrame const* plotting_frame = plugin->renderer().GetPlottingFrame();
  auto const* plotting_frame_as_rotating_pulsating =
      dynamic_cast<
          RotatingPulsatingReferenceFrame<Barycentric,
                                          Navigation> const*>(plotting_frame);
  if (!plotting_frame_as_rotating_pulsating) {
    // We do not draw equipotentials in the current plotting frame.
    return m.Return(0);
  } else {
    // TODO(egg): We may not want to recompute all the time.
    plugin->geometric_potential_plotter().RequestEquipotentials(
        {plotting_frame_as_rotating_pulsating->primaries(),
         plotting_frame_as_rotating_pulsating->secondaries(),
         plugin->CurrentTime()});
    if (equipotentials != nullptr &&
        plotting_frame_as_rotating_pulsating->primaries() ==
            equipotentials->parameters.primaries &&
        plotting_frame_as_rotating_pulsating->secondaries() ==
            equipotentials->parameters.secondaries) {
      return m.Return(equipotentials->lines.size());
    } else {
      return m.Return(0);
    }
  }
}

void __cdecl principia__FreeVesselsAndPartsAndCollectPileUps(
    Plugin* const plugin,
    double const delta_t) {
  journal::Method<journal::FreeVesselsAndPartsAndCollectPileUps> m(
      {plugin, delta_t});
  CHECK_NOTNULL(plugin);
  plugin->FreeVesselsAndPartsAndCollectPileUps(delta_t * Second);
  return m.Return();
}

int __cdecl principia__GetBufferDuration() {
  journal::Method<journal::GetBufferDuration> m;
  return m.Return(FLAGS_logbufsecs);
}

int __cdecl principia__GetBufferedLogging() {
  journal::Method<journal::GetBufferedLogging> m;
  return m.Return(FLAGS_logbuflevel);
}

int __cdecl principia__GetStderrLogging() {
  journal::Method<journal::GetStderrLogging> m;
  return m.Return(FLAGS_stderrthreshold);
}

int __cdecl principia__GetSuppressedLogging() {
  journal::Method<journal::GetSuppressedLogging> m;
  return m.Return(FLAGS_minloglevel);
}

int __cdecl principia__GetVerboseLogging() {
  journal::Method<journal::GetVerboseLogging> m;
  return m.Return(FLAGS_v);
}

void __cdecl principia__GetVersion(
    char const** const build_date,
    char const** const version) {
  journal::Method<journal::GetVersion> m({build_date, version});
  *CHECK_NOTNULL(build_date) = BuildDate;
  *CHECK_NOTNULL(version) = Version;
  return m.Return();
}

bool __cdecl principia__HasEncounteredApocalypse(
    Plugin* const plugin,
    char const** const details) {
  journal::Method<journal::HasEncounteredApocalypse> m({plugin}, {details});
  // Ownership will be transfered to the marshmallow.
  std::string details_string;
  bool const has_encountered_apocalypse =
      CHECK_NOTNULL(plugin)->HasEncounteredApocalypse(&details_string);
  UniqueArray<char> allocated_details(details_string.size() + 1);
  std::memcpy(allocated_details.data.get(),
              details_string.data(),
              details_string.size() + 1);
  *CHECK_NOTNULL(details) = allocated_details.data.release();
  return m.Return(has_encountered_apocalypse);
}

bool __cdecl principia__HasVessel(Plugin* const plugin,
                                  char const* const vessel_guid) {
  journal::Method<journal::HasVessel> m({plugin,  vessel_guid});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(vessel_guid);
  return m.Return(plugin->HasVessel(vessel_guid));
}

// Sets stderr to log INFO, and redirects stderr, which Unity does not log, to
// "<KSP directory>/stderr.log".  This provides an easily accessible file
// containing a sufficiently verbose log of the latest session, instead of
// requiring users to dig in the archive of all past logs at all severities.
// This archive is written to
// "<KSP directory>/glog/Principia/<SEVERITY>.<date>-<time>.<pid>",
// where date and time are in ISO 8601 basic format.
void __cdecl principia__InitGoogleLogging() {
  if (google::IsGoogleLoggingInitialized()) {
    LOG(INFO) << "Google logging was already initialized, no action taken";
  } else {
#ifdef _MSC_VER
    FILE* file;
    freopen_s(&file, "stderr.log", "w", stderr);
#else
    std::freopen("stderr.log", "w", stderr);
#endif
    google::SetLogFilenameExtension(".log");
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

    LOG(ERROR) << "Initialized Google logging for Principia";
    LOG(ERROR) << "Principia version " << Version
               << " built on " << BuildDate
               << " by " << principia::base::CompilerName
               << " version " << principia::base::CompilerVersion
               << " for " << principia::base::OperatingSystem << " "
               << principia::base::Architecture;
    LOG(ERROR) << "Running on " << ProcessorBrandString() << " ("
               << CPUVendorIdentificationString() << ")";
    LOG(ERROR) << "with " << CPUFeatures();
#if OS_WIN
  MODULEINFO module_info;
  memset(&module_info, 0, sizeof(module_info));
  CHECK(GetModuleInformation(GetCurrentProcess(),
                             GetModuleHandle(TEXT("principia")),
                             &module_info,
                             sizeof(module_info)));
  LOG(ERROR) << "Base address is " << module_info.lpBaseOfDll;
#endif
  }
}

void __cdecl principia__InitializeDownsamplingParameters(
    Plugin* const plugin,
    ConfigurationDownsamplingParameters const& downsampling_parameters) {
  journal::Method<journal::InitializeDownsamplingParameters> m(
      {plugin, downsampling_parameters});
  CHECK_NOTNULL(plugin);
  plugin->InitializeDownsamplingParameters(
      MakeDownsamplingParameters(downsampling_parameters));
  return m.Return();
}

void __cdecl principia__InitializeEphemerisParameters(
    Plugin* const plugin,
    ConfigurationAccuracyParameters const& accuracy_parameters,
    ConfigurationFixedStepParameters const& fixed_step_parameters) {
  journal::Method<journal::InitializeEphemerisParameters> m(
      {plugin, accuracy_parameters, fixed_step_parameters});
  CHECK_NOTNULL(plugin);
  plugin->InitializeEphemerisParameters(
      MakeAccuracyParameters(accuracy_parameters),
      MakeFixedStepParameters(fixed_step_parameters));
  return m.Return();
}

void __cdecl principia__InitializeHistoryParameters(
    Plugin* const plugin,
    ConfigurationFixedStepParameters const& fixed_step_parameters) {
  journal::Method<journal::InitializeHistoryParameters> m(
      {plugin, fixed_step_parameters});
  CHECK_NOTNULL(plugin);
  plugin->InitializeHistoryParameters(
      MakeFixedStepParameters(fixed_step_parameters));
  return m.Return();
}

void __cdecl principia__InitializePsychohistoryParameters(
    Plugin* const plugin,
    ConfigurationAdaptiveStepParameters const& parameters) {
  journal::Method<journal::InitializePsychohistoryParameters> m(
      {plugin, parameters});
  CHECK_NOTNULL(plugin);
  plugin->InitializePsychohistoryParameters(
      MakeAdaptiveStepParameters(parameters));
  return m.Return();
}

void __cdecl principia__InsertCelestialAbsoluteCartesian(
    Plugin* const plugin,
    int const celestial_index,
    int const* const parent_index,
    BodyParameters const& body_parameters,
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
       body_parameters,
       x, y, z,
       vx, vy, vz});
  CHECK_NOTNULL(plugin);
  serialization::InitialState::Cartesian::Body initial_state;
  initial_state.set_name(body_parameters.name);
  initial_state.set_x(x);
  initial_state.set_y(y);
  initial_state.set_z(z);
  initial_state.set_vx(vx);
  initial_state.set_vy(vy);
  initial_state.set_vz(vz);
  plugin->InsertCelestialAbsoluteCartesian(
      celestial_index,
      parent_index == nullptr ? std::nullopt
                              : std::make_optional(*parent_index),
      MakeGravityModel(body_parameters),
      initial_state);
  return m.Return();
}

void __cdecl principia__InsertCelestialJacobiKeplerian(
    Plugin* const plugin,
    int const celestial_index,
    int const* const parent_index,
    BodyParameters const& body_parameters,
    KeplerianElements const* const keplerian_elements) {
  journal::Method<journal::InsertCelestialJacobiKeplerian> m(
      {plugin,
       celestial_index,
       parent_index,
       body_parameters,
       keplerian_elements});
  CHECK_NOTNULL(plugin);
  serialization::InitialState::Keplerian::Body initial_state;
  initial_state.set_name(body_parameters.name);
  if (keplerian_elements != nullptr) {
    serialization::InitialState::Keplerian::Body::Elements* elements =
        initial_state.mutable_elements();
    elements->set_eccentricity(keplerian_elements->eccentricity);
    if (!std::isnan(keplerian_elements->semimajor_axis)) {
      elements->set_semimajor_axis(
          DebugString(keplerian_elements->semimajor_axis * Metre));
    }
    if (!std::isnan(keplerian_elements->mean_motion)) {
      elements->set_mean_motion(
          DebugString(keplerian_elements->mean_motion * Radian / Second));
    }
    elements->set_inclination(
        DebugString(keplerian_elements->inclination_in_degrees * Degree));
    elements->set_longitude_of_ascending_node(DebugString(
        keplerian_elements->longitude_of_ascending_node_in_degrees * Degree));
    elements->set_argument_of_periapsis(DebugString(
        keplerian_elements->argument_of_periapsis_in_degrees * Degree));
    elements->set_mean_anomaly(
        DebugString(keplerian_elements->mean_anomaly * Radian));
  }
  plugin->InsertCelestialJacobiKeplerian(
      celestial_index,
      parent_index == nullptr ? std::nullopt
                              : std::make_optional(*parent_index),
      MakeGravityModel(body_parameters),
      initial_state);
  return m.Return();
}

// Calls `plugin->InsertOrKeepVessel` with the arguments given.
// `plugin` must not be null.  No transfer of ownership.
void __cdecl principia__InsertOrKeepVessel(Plugin* const plugin,
                                           char const* const vessel_guid,
                                           char const* const vessel_name,
                                           int const parent_index,
                                           bool const loaded,
                                           bool* inserted) {
  journal::Method<journal::InsertOrKeepVessel> m(
      {plugin, vessel_guid, vessel_name, parent_index, loaded}, {inserted});
  CHECK_NOTNULL(plugin);
  CHECK(vessel_guid != nullptr)
      << "name: "
      << (vessel_name == nullptr ? "null" : std::string_view(vessel_name));
  CHECK(vessel_name != nullptr)
      << "guid: "
      << (vessel_guid == nullptr ? "null" : std::string_view(vessel_guid));
  plugin->InsertOrKeepVessel(vessel_guid,
                             vessel_name,
                             parent_index,
                             loaded,
                             *inserted);
  return m.Return();
}

void __cdecl principia__InsertOrKeepLoadedPart(
    Plugin* const plugin,
    PartId const part_id,
    char const* const name,
    double const mass_in_tonnes,
    XYZ const centre_of_mass,
    XYZ const moments_of_inertia_in_tonnes,
    WXYZ const principal_axes_rotation,
    bool const is_solid_rocket_motor,
    char const* const vessel_guid,
    int const main_body_index,
    QP const main_body_world_degrees_of_freedom,
    QP const part_world_degrees_of_freedom,
    WXYZ const part_rotation,
    XYZ const part_angular_velocity,
    double const delta_t) {
  journal::Method<journal::InsertOrKeepLoadedPart> m(
      {plugin,
       part_id,
       name,
       mass_in_tonnes,
       centre_of_mass,
       moments_of_inertia_in_tonnes,
       principal_axes_rotation,
       is_solid_rocket_motor,
       vessel_guid,
       main_body_index,
       main_body_world_degrees_of_freedom,
       part_world_degrees_of_freedom,
       part_rotation,
       part_angular_velocity,
       delta_t});
  CHECK_NOTNULL(plugin);

  // We build the inertia tensor in the principal axes and then transform it to
  // RigidPart.
  using PartPrincipalAxes = Frame<serialization::Frame::PhysicsTag,
                                  Arbitrary,
                                  Handedness::Left,
                                  serialization::Frame::PRINCIPAL_AXES>;

  static constexpr MomentOfInertia zero;

  auto const moments_of_inertia =
      FromXYZ<R3Element<MomentOfInertia>>({moments_of_inertia_in_tonnes.x,
                                           moments_of_inertia_in_tonnes.y,
                                           moments_of_inertia_in_tonnes.z});
  InertiaTensor<PartPrincipalAxes> const inertia_tensor_in_princial_axes(
      R3x3Matrix<MomentOfInertia>({moments_of_inertia.x, zero, zero},
                                  {zero, moments_of_inertia.y, zero},
                                  {zero, zero, moments_of_inertia.z}));

  Rotation<PartPrincipalAxes, RigidPart> const principal_axes_to_rigid_part(
      FromWXYZ(principal_axes_rotation));
  InertiaTensor<RigidPart> const inertia_tensor_in_rigid_part =
      principal_axes_to_rigid_part(inertia_tensor_in_princial_axes);

  VLOG(1) << "InsertOrKeepLoadedPart: " << name << " " << part_id << " "
          << moments_of_inertia << " " << FromWXYZ(principal_axes_rotation);

  plugin->InsertOrKeepLoadedPart(
      part_id,
      name,
      mass_in_tonnes * Tonne,
      FromXYZ<Position<EccentricPart>>(centre_of_mass),
      inertia_tensor_in_rigid_part,
      is_solid_rocket_motor,
      vessel_guid,
      main_body_index,
      FromQP<DegreesOfFreedom<World>>(main_body_world_degrees_of_freedom),
      MakePartRigidMotion(
          part_world_degrees_of_freedom, part_rotation, part_angular_velocity),
      delta_t * Second);
  return m.Return();
}

// Calls `plugin->SetVesselStateOffset` with the arguments given.
// `plugin` must not be null.  No transfer of ownership.
void __cdecl principia__InsertUnloadedPart(Plugin* const plugin,
                                           PartId const part_id,
                                           char const* const name,
                                           char const* const vessel_guid,
                                           QP const from_parent) {
  journal::Method<journal::InsertUnloadedPart> m(
      {plugin, part_id, name, vessel_guid, from_parent});
  CHECK_NOTNULL(plugin);
  plugin->InsertUnloadedPart(
      part_id,
      name,
      vessel_guid,
      FromQP<RelativeDegreesOfFreedom<AliceSun>>(from_parent));
  return m.Return();
}

// Exports `LOG(SEVERITY) << text` for fast logging from the C# adapter.
// This will always evaluate its argument even if the corresponding log severity
// is disabled, so it is less efficient than LOG(SEVERITY).
void __cdecl principia__LogError(char const* const file,
                                 int const line,
                                 char const* const text) {
  journal::Method<journal::LogError> m({file, line, text});
  google::LogMessage(file, line, google::ERROR).stream() << text;
  return m.Return();
}

void __cdecl principia__LogFatal(char const* const file,
                                 int const line,
                                 char const* const text) {
  journal::Method<journal::LogFatal> m({file, line, text});
  google::LogMessageFatal(file, line).stream() << text;
  return m.Return();
}

void __cdecl principia__LogInfo(char const* const file,
                                int const line,
                                char const* const text) {
  journal::Method<journal::LogInfo> m({file, line, text});
  google::LogMessage(file, line).stream() << text;
  return m.Return();
}

void __cdecl principia__LogWarning(char const* const file,
                                   int const line,
                                   char const* const text) {
  journal::Method<journal::LogWarning> m({file, line, text});
  google::LogMessage(file, line, google::WARNING).stream() << text;
  return m.Return();
}

WXYZ __cdecl principia__NavballOrientation(
    Plugin const* const plugin,
    XYZ const sun_world_position,
    XYZ const ship_world_position) {
  journal::Method<journal::NavballOrientation> m({plugin,
                                                  sun_world_position,
                                                  ship_world_position});
  CHECK_NOTNULL(plugin);
  auto const frame_field = plugin->NavballFrameField(
      FromXYZ<Position<World>>(sun_world_position));
  return m.Return(ToWXYZ(
      frame_field->FromThisFrame(
          FromXYZ<Position<World>>(ship_world_position)).quaternion()));
}

// Returns a pointer to a plugin constructed with the arguments given.
// The caller takes ownership of the result.
Plugin* __cdecl principia__NewPlugin(
    char const* const game_epoch,
    char const* const solar_system_epoch,
    double const planetarium_rotation_in_degrees) {
  journal::Method<journal::NewPlugin> m({game_epoch,
                                         solar_system_epoch,
                                         planetarium_rotation_in_degrees});
  LOG(INFO) << "Constructing Principia plugin";
  not_null<std::unique_ptr<Plugin>> result =
      make_not_null_unique<Plugin>(game_epoch,
                                   solar_system_epoch,
                                   planetarium_rotation_in_degrees * Degree);
  LOG(INFO) << "Plugin constructed";
  return m.Return(result.release());
}

void __cdecl principia__PrepareToReportCollisions(Plugin* const plugin) {
  journal::Method<journal::PrepareToReportCollisions> m({plugin});
  CHECK_NOTNULL(plugin)->PrepareToReportCollisions();
  return m.Return();
}

void __cdecl principia__ReportGroundCollision(Plugin const* const plugin,
                                              uint32_t const part_id) {
  journal::Method<journal::ReportGroundCollision> m({plugin, part_id});
  CHECK_NOTNULL(plugin)->ReportGroundCollision(part_id);
  return m.Return();
}

void __cdecl principia__ReportPartCollision(Plugin const* const plugin,
                                            PartId const part1_id,
                                            PartId const part2_id) {
  journal::Method<journal::ReportPartCollision> m({plugin, part1_id, part2_id});
  CHECK_NOTNULL(plugin)->ReportPartCollision(part1_id, part2_id);
  return m.Return();
}

// Says hello, convenient for checking that calls to the DLL work.
char const* __cdecl principia__SayHello() {
  journal::Method<journal::SayHello> m;
  return m.Return("Hello from native C++!");
}

// For checking that interchange messages work.
Status* __cdecl principia__SayNotFound() {
  journal::Method<journal::SayNotFound> m;
  return m.Return(ToNewStatus(
      absl::NotFoundError("Not found from native C++!")));
}

#define PRINCIPIA_VERIFY_SERIALIZATION 0

#if PRINCIPIA_VERIFY_SERIALIZATION
static PushDeserializer* verification_deserializer = nullptr;
static Plugin const* verification_plugin = nullptr;
#endif

// `plugin` must not be null.  The caller takes ownership of the result, except
// when it is null (at the end of the stream).  No transfer of ownership of
// `*plugin`.  `*serializer` must be null on the first call and must be passed
// unchanged to the successive calls; its ownership is not transferred.
char const* __cdecl principia__SerializePlugin(
    Plugin const* const plugin,
    PullSerializer** const serializer,
    char const* const compressor,
    char const* const encoder) {
  journal::Method<journal::SerializePlugin> m({plugin,
                                               serializer,
                                               compressor,
                                               encoder},
                                              {serializer});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(serializer);

  // Create and start a serializer if the caller didn't provide one.
  if (*serializer == nullptr) {
    LOG(INFO) << "Begin plugin serialization";
    *serializer = new PullSerializer(chunk_size,
                                     number_of_chunks,
                                     NewCompressor(compressor));
    not_null<serialization::Plugin*> const message =
        Arena::CreateMessage<serialization::Plugin>(arena);
    plugin->WriteToMessage(message);
    (*serializer)->Start(message);
  }

  // Pull a chunk.
  Array<std::uint8_t> bytes;
  bytes = (*serializer)->Pull();

  // If this is the end of the serialization, delete the serializer and return a
  // nullptr.
  if (bytes.size == 0) {
#if PRINCIPIA_VERIFY_SERIALIZATION
    principia__DeserializePlugin("",
                                 &verification_deserializer,
                                 &verification_plugin,
                                 compressor,
                                 encoder);
    CHECK_NOTNULL(verification_plugin);
    LOG(INFO) << "Deleting verification plugin";
    delete verification_plugin;
    verification_plugin = nullptr;
    LOG(INFO) << "Verification plugin deleted";
#endif
    LOG(INFO) << "End plugin serialization";
    TakeOwnership(serializer);
    arena->Reset();
    return m.Return(nullptr);
  }

  // Encode and return to the client.
  auto hexadecimal = NewEncoder(encoder)->Encode(bytes);
#if PRINCIPIA_VERIFY_SERIALIZATION
  principia__DeserializePlugin(hexadecimal.data.get(),
                               &verification_deserializer,
                               &verification_plugin,
                               compressor,
                               encoder);
#endif
  return m.Return(hexadecimal.data.release());
}

// Sets the maximum number of seconds which logs may be buffered for.
void __cdecl principia__SetBufferDuration(int const seconds) {
  journal::Method<journal::SetBufferDuration> m({seconds});
  FLAGS_logbufsecs = seconds;
  return m.Return();
}

// Log messages at a level `<= max_severity` are buffered.
// Log messages at a higher level are flushed immediately.
void __cdecl principia__SetBufferedLogging(int const max_severity) {
  journal::Method<journal::SetBufferedLogging> m({max_severity});
  FLAGS_logbuflevel = max_severity;
  return m.Return();
}

void __cdecl principia__SetFlag(char const* const name,
                                char const* const value) {
  journal::Method<journal::SetFlag> m({name, value});
  Flags::Set(name, value);
  return m.Return();
}

void __cdecl principia__SetMainBody(Plugin* const plugin, int const index) {
  journal::Method<journal::SetMainBody> m({plugin, index});
  CHECK_NOTNULL(plugin);
  plugin->SetMainBody(index);
  return m.Return();
}

// Make it so that all log messages of at least `min_severity` are logged to
// stderr (in addition to logging to the usual log file(s)).
void __cdecl principia__SetStderrLogging(int const min_severity) {
  journal::Method<journal::SetStderrLogging> m({min_severity});
  // NOTE(egg): We could use `FLAGS_stderrthreshold` instead, the difference
  // seems to be a mutex.
  google::SetStderrLogging(min_severity);
  return m.Return();
}

// Log suppression level: messages logged at a lower level than this are
// suppressed.
void __cdecl principia__SetSuppressedLogging(int const min_severity) {
  journal::Method<journal::SetSuppressedLogging> m({min_severity});
  FLAGS_minloglevel = min_severity;
  return m.Return();
}

// Show all VLOG(m) messages for `m <= level`.
void __cdecl principia__SetVerboseLogging(int const level) {
  journal::Method<journal::SetVerboseLogging> m({level});
  FLAGS_v = level;
  return m.Return();
}

void __cdecl principia__SetWorldRotationalReferenceFrame(Plugin* const plugin,
                                                         int const index) {
  journal::Method<journal::SetWorldRotationalReferenceFrame> m({plugin, index});
  CHECK_NOTNULL(plugin);
  plugin->SetWorldRotationalReferenceFrame(index);
  return m.Return();
}

XYZ __cdecl principia__UnmanageableVesselVelocity(Plugin const* const plugin,
                                                  QP const degrees_of_freedom,
                                                  int const parent_index) {
  return ToXYZ(CHECK_NOTNULL(plugin)->UnmanageableVesselVelocity(
      FromQP<RelativeDegreesOfFreedom<AliceSun>>(degrees_of_freedom),
      parent_index));
}

// Calls `plugin->UpdateCelestialHierarchy` with the arguments given.
// `plugin` must not be null.  No transfer of ownership.
void __cdecl principia__UpdateCelestialHierarchy(Plugin const* const plugin,
                                                 int const celestial_index,
                                                 int const parent_index) {
  journal::Method<journal::UpdateCelestialHierarchy> m({plugin,
                                                        celestial_index,
                                                        parent_index});
  CHECK_NOTNULL(plugin);
  plugin->UpdateCelestialHierarchy(celestial_index, parent_index);
  return m.Return();
}

void __cdecl principia__UpdatePrediction(
    Plugin const* const plugin,
    char const* const* const vessel_guids) {
  journal::Method<journal::UpdatePrediction> m({plugin, vessel_guids});
  CHECK_NOTNULL(plugin);
  std::vector<GUID> guids;
  for (char const* const* c = vessel_guids;
       *c != nullptr;
       ++c) {
    guids.push_back(std::string(*c));
  }
  plugin->UpdatePrediction(guids);
  return m.Return();
}

}  // namespace interface
}  // namespace principia
