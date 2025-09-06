#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "astronomy/solar_system_fingerprints.hpp"
#include "astronomy/stabilize_ksp.hpp"
#include "astronomy/time_scales.hpp"
#include "base/file.hpp"
#include "base/fingerprint2011.hpp"
#include "base/flags.hpp"
#include "base/hexadecimal.hpp"
#include "base/map_util.hpp"
#include "base/serialization.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/identity.hpp"
#include "geometry/permutation.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/sign.hpp"
#include "geometry/space_transformations.hpp"
#include "glog/logging.h"
#include "glog/stl_logging.h"
#include "ksp_plugin/equator_relevance_threshold.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/part_subsets.hpp"  // 🧙 For Subset<Part>.
#include "physics/apsides.hpp"
#include "physics/barycentric_rotating_reference_frame.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/body_surface_frame_field.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/parser.hpp"

namespace principia {
namespace ksp_plugin {
namespace _plugin {
namespace internal {

using ::operator<<;
using namespace principia::astronomy::_solar_system_fingerprints;
using namespace principia::astronomy::_stabilize_ksp;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_file;
using namespace principia::base::_fingerprint2011;
using namespace principia::base::_flags;
using namespace principia::base::_hexadecimal;
using namespace principia::base::_map_util;
using namespace principia::base::_serialization;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_identity;
using namespace principia::geometry::_permutation;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_space_transformations;
using namespace principia::ksp_plugin::_equator_relevance_threshold;
using namespace principia::ksp_plugin::_integrators;
using namespace principia::ksp_plugin::_part;
using namespace principia::physics::_apsides;
using namespace principia::physics::_barycentric_rotating_reference_frame;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_body_surface_frame_field;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_parser;

Length const& MaxCollisionError() {
  static Length const max_collision_error = []() {
    std::string_view name = "max_collision_error";
    if (Flags::IsPresent(name)) {
      auto const values = Flags::Values(name);
      CHECK_EQ(values.size(), 1);
      return ParseQuantity<Length>(
          *Flags::Values("max_collision_error").begin());
    } else {
      return 10 * Metre;
    }
  }();
  return max_collision_error;
}

// Keep this consistent with `prediction_steps_` in `main_window.cs`.
constexpr std::int64_t max_steps_in_prediction = 1 << 24;

Plugin::Plugin(std::string const& game_epoch,
               std::string const& solar_system_epoch,
               Angle const& planetarium_rotation)
    : history_downsampling_parameters_(DefaultDownsamplingParameters()),
      history_fixed_step_parameters_(DefaultHistoryParameters()),
      psychohistory_parameters_(DefaultPsychohistoryParameters()),
      vessel_thread_pool_(
          /*pool_size=*/2 * std::thread::hardware_concurrency()),
      planetarium_rotation_(planetarium_rotation),
      game_epoch_(ParseTT(game_epoch)),
      current_time_(ParseTT(solar_system_epoch)) {
  gravity_model_.set_plugin_frame(serialization::Frame::BARYCENTRIC);
  initial_state_.set_epoch(solar_system_epoch);
  initial_state_.set_plugin_frame(serialization::Frame::BARYCENTRIC);
}

Plugin::~Plugin() {
  // We must manually destroy the vessels, triggering the destruction of the
  // parts, which have callbacks to remove themselves from `part_id_to_vessel_`,
  // which must therefore still exist.  This also causes the parts to be
  // destroyed, and therefore to destroy the pile-ups, which want to remove
  // themselves from `pile_up_`, which also exists.
  vessels_.clear();
}

void Plugin::InsertCelestialAbsoluteCartesian(
    Index const celestial_index,
    std::optional<Index> const& parent_index,
    serialization::GravityModel::Body const& gravity_model,
    serialization::InitialState::Cartesian::Body const& initial_state) {
  CHECK_EQ(gravity_model.name(), initial_state.name());
  InitializeIndices(gravity_model.name(), celestial_index, parent_index);
  *gravity_model_.add_body() = gravity_model;
  CHECK(!initial_state_.has_keplerian()) << initial_state_.DebugString();
  *initial_state_.mutable_cartesian()->add_body() = initial_state;
}

void Plugin::InsertCelestialJacobiKeplerian(
    Index const celestial_index,
    std::optional<Index> const& parent_index,
    serialization::GravityModel::Body const& gravity_model,
    serialization::InitialState::Keplerian::Body const& initial_state) {
  CHECK_EQ(gravity_model.name(), initial_state.name());
  InitializeIndices(gravity_model.name(), celestial_index, parent_index);
  *gravity_model_.add_body() = gravity_model;
  CHECK(!initial_state_.has_cartesian()) << initial_state_.DebugString();
  serialization::InitialState::Keplerian::Body* const body =
      initial_state_.mutable_keplerian()->add_body();
  *body = initial_state;
  if (parent_index) {
    body->set_parent(FindOrDie(index_to_name_, *parent_index));
  }
}

void Plugin::InitializeDownsamplingParameters(
    DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters const&
        downsampling_parameters) {
  CHECK(initializing_);
  history_downsampling_parameters_ = downsampling_parameters;
}

void Plugin::InitializeEphemerisParameters(
    Ephemeris<Barycentric>::AccuracyParameters const& accuracy_parameters,
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters) {
  CHECK(initializing_);
  ephemeris_accuracy_parameters_ = accuracy_parameters;
  ephemeris_fixed_step_parameters_ = fixed_step_parameters;
}

void Plugin::InitializeHistoryParameters(
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters) {
  CHECK(initializing_);
  history_fixed_step_parameters_ = fixed_step_parameters;
}

void Plugin::InitializePsychohistoryParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const& parameters) {
  CHECK(initializing_);
  psychohistory_parameters_ = parameters;
}

void Plugin::EndInitialization() {
  CHECK(initializing_);
  SolarSystem<Barycentric> solar_system(gravity_model_, initial_state_);

  // Check if this is the stock KSP system in which case it needs to be
  // stabilized.
  system_fingerprint_ = solar_system.Fingerprint();
  LOG(INFO) << "System fingerprint is 0x" << std::hex << std::uppercase
            << system_fingerprint_;

  bool is_well_known = false;
  for (auto const ksp_version : {KSP122, KSP191}) {
    if (system_fingerprint_ == KSPStockSystemFingerprints[ksp_version]) {
      LOG(WARNING) << "This appears to be the dreaded KSP stock system!";
      StabilizeKSP(solar_system);
      system_fingerprint_ = solar_system.Fingerprint();
      LOG(INFO) << "System fingerprint after stabilization is 0x" << std::hex
                << std::uppercase << system_fingerprint_;
      CHECK_EQ(KSPStabilizedSystemFingerprints[ksp_version],
               system_fingerprint_)
          << "Attempt at stabilizing the KSP system failed!\n"
          << gravity_model_.DebugString() << "\n"
          << initial_state_.DebugString();
      LOG(INFO) << "This is the stabilized KSP system, all hail retrobop!";
      is_well_known = true;
      break;
    } else if (system_fingerprint_ ==
                KSPStabilizedSystemFingerprints[ksp_version]) {
      LOG(INFO) << "This is the stabilized KSP system, and we didn't have to "
                << "stabilize it ourselves.  All hail retrobop anyway!";
      is_well_known = true;
      break;
    }
  }
  if (!is_well_known) {
    LOG(WARNING) << "This is an unknown system, we don't know anything about "
                 << "its stability:\n"
                 << gravity_model_.DebugString() << "\n"
                 << initial_state_.DebugString();
  }

  // Construct the ephemeris.
  ephemeris_ =
      solar_system.MakeEphemeris(ephemeris_accuracy_parameters_.value_or(
                                     DefaultEphemerisAccuracyParameters()),
                                 ephemeris_fixed_step_parameters_.value_or(
                                     DefaultEphemerisFixedStepParameters()));

  // Construct the celestials using the bodies from the ephemeris.
  for (std::string const& name : solar_system.names()) {
    auto const rotating_body = solar_system.rotating_body(*ephemeris_, name);
    Index const celestial_index = FindOrDie(name_to_index_, name);
    auto const [it, inserted] = celestials_.emplace(
        celestial_index, std::make_unique<Celestial>(rotating_body));
    CHECK(inserted) << "Body already exists at index " << celestial_index;
    it->second->set_trajectory(ephemeris_->trajectory(rotating_body));
  }

  // Establish the parent relationships between the celestials.
  for (auto const& [celestial_index, celestial] : celestials_) {
    auto const& parent_index = FindOrDie(parents_, celestial_index);
    if (parent_index) {
      not_null<Celestial const*> parent =
          FindOrDie(celestials_, *parent_index).get();
      celestial->set_parent(parent);
    } else {
      CHECK(sun_ == nullptr);
      sun_ = celestial.get();
    }
  }
  CHECK_NOTNULL(sun_);
  main_body_ = sun_->body();

  UpdatePlanetariumRotation();

  // This would use NewBodyCentredNonRotatingNavigationFrame, but we don't have
  // the sun's index at hand.
  // TODO(egg): maybe these functions should take `Celestial*`s, and we should
  // then export `FindOrDie(celestials_, _)`.
  renderer_ = std::make_unique<Renderer>(
      sun_,
      make_not_null_unique<
          BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation>>(
              ephemeris_.get(),
              sun_->body()));

  geometric_potential_plotter_.emplace(ephemeris_.get());

  // Log the serialized ephemeris.
  serialization::Ephemeris ephemeris_message;
  ephemeris_->WriteToMessage(&ephemeris_message);
  HexadecimalEncoder</*null_terminated=*/true> encoder;
  auto const hex = encoder.Encode(SerializeAsBytes(ephemeris_message).get());
  // Begin and end markers to make sure the hex did not get clipped (this might
  // happen if the message is very big).
  LOG(INFO) << "Ephemeris at initialization:\nbegin\n"
            << hex.data.get() << "\nend";

  initializing_.Flop();
}

bool Plugin::HasEncounteredApocalypse(std::string* const details) const {
  CHECK_NOTNULL(details);
  auto const status = ephemeris_->last_severe_integration_status();
  if (absl::IsInvalidArgument(status)) {
    *details = status.message();
    return true;
  } else {
    return false;
  }
}

void Plugin::UpdateCelestialHierarchy(Index const celestial_index,
                                      Index const parent_index) const {
  CHECK(!initializing_);
  FindOrDie(celestials_, celestial_index)->set_parent(
      FindOrDie(celestials_, parent_index).get());
}

void Plugin::SetMainBody(Index const index) {
  main_body_ = FindOrDie(celestials_, index)->body();
  LOG_IF(FATAL, main_body_ == nullptr) << index;
  UpdatePlanetariumRotation();
}

Rotation<BodyWorld, World> Plugin::CelestialRotation(
    Index const index) const {
  // `BodyWorld` with its y and z axes swapped (so that z is the polar axis).
  using BodyFixed = Frame<struct BodyFixedTag>;
  Permutation<BodyWorld, BodyFixed> const body_mirror(OddPermutation::XZY);

  auto const& body = *FindOrDie(celestials_, index)->body();

  OrthogonalMap<BodyWorld, World> const result =
      OrthogonalMap<WorldSun, World>::Identity() *
      sun_looking_glass.Inverse().Forget<OrthogonalMap>() *
      (PlanetariumRotation() * body.FromSurfaceFrame<BodyFixed>(current_time_))
          .Forget<OrthogonalMap>() *
      body_mirror.Forget<OrthogonalMap>();
  return result.AsRotation();
}

Rotation<CelestialSphere, World> Plugin::CelestialSphereRotation()
    const {
  Permutation<CelestialSphere, Barycentric> const celestial_mirror(
      OddPermutation::XZY);
  auto const result = OrthogonalMap<WorldSun, World>::Identity() *
                      sun_looking_glass.Inverse().Forget<OrthogonalMap>() *
                      PlanetariumRotation().Forget<OrthogonalMap>() *
                      celestial_mirror.Forget<OrthogonalMap>();
  return result.AsRotation();
}

Angle Plugin::CelestialInitialRotation(Index const celestial_index) const {
  auto const& body = *FindOrDie(celestials_, celestial_index)->body();
  // Offset by π/2 since `AngleAt` is with respect to the y axis of the
  // celestial frame of the body, but KSP counts from the x axis.
  return body.AngleAt(game_epoch_) + π / 2 * Radian;
}

Time Plugin::CelestialRotationPeriod(Index const celestial_index) const {
  auto const& body = *FindOrDie(celestials_, celestial_index)->body();
  // The result will be negative if the pole is the negative pole
  // (e.g. for Venus).  This is the convention KSP uses for retrograde rotation.
  return 2 * π * Radian / body.angular_frequency();
}

void Plugin::ClearWorldRotationalReferenceFrame() {
  angular_velocity_of_world_ = Barycentric::nonrotating;
}

void Plugin::SetWorldRotationalReferenceFrame(Index const celestial_index) {
  SetMainBody(celestial_index);
  angular_velocity_of_world_ = main_body_->angular_velocity();
}

Index Plugin::CelestialIndexOfBody(MassiveBody const& body) const {
  return FindOrDie(name_to_index_, body.name());
}

void Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                std::string const& vessel_name,
                                Index const parent_index,
                                bool const loaded,
                                bool& inserted) {
  CHECK(!initializing_);
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  auto vit = vessels_.find(vessel_guid);
  if (vit == vessels_.end()) {
    // Restore the zombie parameters if we have some, otherwise use the default.
    auto prediction_parameters = DefaultPredictionParameters();
    if (auto const pit =
            zombie_prediction_adaptive_step_parameters_.find(vessel_guid);
        pit != zombie_prediction_adaptive_step_parameters_.end()) {
      prediction_parameters = pit->second;
      zombie_prediction_adaptive_step_parameters_.erase(pit);
    }
    std::tie(vit, inserted) =
        vessels_.emplace(
            vessel_guid,
            make_not_null_unique<Vessel>(vessel_guid,
                                         vessel_name,
                                         parent,
                                         ephemeris_.get(),
                                         prediction_parameters,
                                         history_downsampling_parameters_));
  } else {
    inserted = false;
  }
  not_null<Vessel*> const vessel = vit->second.get();
  if (vessel->name() != vessel_name) {
    vessel->set_name(vessel_name);
  }
  kept_vessels_.emplace(vessel);
  if (vessel->parent() != parent) {
    vessel->set_parent(parent);
  }
  if (loaded) {
    loaded_vessels_.insert(vessel);
  }
  LOG_IF(INFO, inserted) << "Inserted " << (loaded ? "loaded" : "unloaded")
                         << " vessel " << vessel->ShortDebugString();
}

void Plugin::InsertUnloadedPart(
    PartId part_id,
    std::string const& name,
    GUID const& vessel_guid,
    RelativeDegreesOfFreedom<AliceSun> const& from_parent) {
  not_null<Vessel*> const vessel = FindOrDie(vessels_, vessel_guid).get();
  ephemeris_->Prolong(current_time_).IgnoreError();

  RelativeDegreesOfFreedom<Barycentric> const relative =
      PlanetariumRotation().Inverse()(from_parent);
  DegreesOfFreedom<Barycentric> const degrees_of_freedom =
      vessel->parent()->current_degrees_of_freedom(current_time_) + relative;

  AddPart(vessel, part_id, name, degrees_of_freedom);
  // NOTE(egg): we do not keep the part; it may disappear just as we load, if
  // it happens to be a part with no physical significance (rb == null).
}

void Plugin::InsertOrKeepLoadedPart(
    PartId const part_id,
    std::string const& name,
    Mass const& mass,
    Position<EccentricPart> const& centre_of_mass,
    InertiaTensor<RigidPart> const& inertia_tensor,
    bool const is_solid_rocket_motor,
    GUID const& vessel_guid,
    Index const main_body_index,
    DegreesOfFreedom<World> const& main_body_degrees_of_freedom,
    RigidMotion<EccentricPart, World> const& part_rigid_motion,
    Time const& Δt) {
  not_null<Vessel*> const vessel = FindOrDie(vessels_, vessel_guid).get();
  CHECK(is_loaded(vessel));

  Instant const previous_time = current_time_ - Δt;
  OrthogonalMap<Barycentric, Barycentric> const Δplanetarium_rotation =
      Exp(Δt * angular_velocity_of_world_).Forget<OrthogonalMap>();
  // TODO(egg): Can we use `BarycentricToWorld` here?
  BodyCentredNonRotatingReferenceFrame<Barycentric, MainBodyCentred> const
      main_body_frame{ephemeris_.get(),
                      FindOrDie(celestials_, main_body_index)->body()};
  RigidMotion<World, MainBodyCentred> const world_to_main_body_centred{
      RigidTransformation<World, MainBodyCentred>{
          main_body_degrees_of_freedom.position(),
          MainBodyCentred::origin,
          main_body_frame.ToThisFrameAtTime(previous_time).orthogonal_map() *
              Δplanetarium_rotation.Inverse() *
              renderer_->WorldToBarycentric(PlanetariumRotation())},
      (renderer_->BarycentricToWorld(PlanetariumRotation()) *
          Δplanetarium_rotation)(-angular_velocity_of_world_),
      main_body_degrees_of_freedom.velocity()};
  RigidMotion<World, Barycentric> const world_to_barycentric_motion =
      main_body_frame.FromThisFrameAtTime(previous_time) *
      world_to_main_body_centred;

  if (auto const it = part_id_to_vessel_.find(part_id);
      it != part_id_to_vessel_.end()) {
    not_null<Vessel*>& associated_vessel = it->second;
    not_null<Vessel*> const current_vessel = associated_vessel;
    if (vessel == current_vessel) {
    } else {
      associated_vessel = vessel;
      vessel->AddPart(current_vessel->ExtractPart(part_id));
    }
  } else {
    AddPart(vessel,
            part_id,
            name,
            mass,
            centre_of_mass,
            inertia_tensor,
            world_to_barycentric_motion * part_rigid_motion);
  }
  vessel->KeepPart(part_id);
  not_null<Part*> part = vessel->part(part_id);
  part->make_truthful();
  part->set_mass(mass);
  part->set_centre_of_mass(centre_of_mass);
  part->set_is_solid_rocket_motor(is_solid_rocket_motor);
  part->set_inertia_tensor(inertia_tensor);
}

void Plugin::ApplyPartIntrinsicForce(PartId const part_id,
                                     Vector<Force, World> const& force) const {
  CHECK(!initializing_);
  not_null<Vessel*> const vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  vessel->part(part_id)->apply_intrinsic_force(
      renderer_->WorldToBarycentric(PlanetariumRotation())(force));
}

void Plugin::ApplyPartIntrinsicForceAtPosition(
    PartId const part_id,
    Vector<Force, World> const& force,
    Displacement<World> const& lever_arm) const {
  CHECK(!initializing_);
  not_null<Vessel*> const vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  OrthogonalMap<World, Barycentric> const world_to_barycentric =
      renderer_->WorldToBarycentric(PlanetariumRotation());
  vessel->part(part_id)->ApplyIntrinsicForceWithLeverArm(
      world_to_barycentric(force),
      world_to_barycentric(lever_arm));
}

void Plugin::ApplyPartIntrinsicTorque(
    PartId const part_id,
    Bivector<Torque, World> const& torque) const {
  CHECK(!initializing_);
  not_null<Vessel*> const vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  vessel->part(part_id)->apply_intrinsic_torque(
      renderer_->WorldToBarycentric(PlanetariumRotation())(torque));
}

bool Plugin::PartIsTruthful(PartId const part_id) const {
  if (auto const it = part_id_to_vessel_.find(part_id);
      it == part_id_to_vessel_.end()) {
    return false;
  } else {
    return it->second->part(part_id)->truthful();
  }
}

void Plugin::PrepareToReportCollisions() {
  for (auto const& [_, vessel] : vessels_) {
    // NOTE(egg): The lifetime requirement on the second argument of
    // `MakeSingleton` (which forwards to the argument of the constructor of
    // `Subset<Part>::Properties`) is that `part` outlives the constructed
    // `Properties`; since these are owned by `part`, this is true.
    vessel->ForAllParts(
        [](Part& part) { Subset<Part>::MakeSingleton(part, &part); });
  }
}

void Plugin::ReportGroundCollision(PartId const part) const {
  Vessel const& v = *FindOrDie(part_id_to_vessel_, part);
  Part& p = *v.part(part);
  LOG(INFO) << "Collision between " << p.ShortDebugString()
            << " and the ground.";
  Subset<Part>::Find(p).mutable_properties().Ground();
}

void Plugin::ReportPartCollision(PartId const part1, PartId const part2) const {
  Vessel const& v1 = *FindOrDie(part_id_to_vessel_, part1);
  Vessel const& v2 = *FindOrDie(part_id_to_vessel_, part2);
  Part& p1 = *v1.part(part1);
  Part& p2 = *v2.part(part2);
  LOG(INFO) << "Collision between " << p1.ShortDebugString() << " and "
            << p2.ShortDebugString();
  CHECK(Contains(kept_vessels_, &v1)) << v1.ShortDebugString()
                                      << " will vanish";
  CHECK(Contains(kept_vessels_, &v2)) << v2.ShortDebugString()
                                      << " will vanish";
  CHECK(v1.WillKeepPart(part1)) << p1.ShortDebugString() << " will vanish";
  CHECK(v2.WillKeepPart(part2)) << p2.ShortDebugString() << " will vanish";
  Subset<Part>::Unite(Subset<Part>::Find(p1), Subset<Part>::Find(p2));
}

void Plugin::FreeVesselsAndPartsAndCollectPileUps(Time const& Δt) {
  CHECK(!initializing_);

  // Remove the vessels that we don't want to keep.  Vessels that are not kept
  // have had no reported collisions, so their part subsets do not intersect
  // with the subsets in kept vessels, and none of the part subsets that remain
  // contain deleted parts.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    not_null<Vessel*> const vessel = it->second.get();
    Instant const vessel_time =
        is_loaded(vessel) ? current_time_ - Δt : current_time_;
    if (kept_vessels_.erase(vessel) > 0) {
      vessel->CreateTrajectoryIfNeeded(vessel_time);
      ++it;
    } else {
      loaded_vessels_.erase(vessel);
      LOG(INFO) << "Removing vessel " << vessel->ShortDebugString();
      renderer_->ClearTargetVesselIf(vessel);
      zombie_prediction_adaptive_step_parameters_.insert_or_assign(
          vessel->guid(), vessel->prediction_adaptive_step_parameters());
      it = vessels_.erase(it);
    }
  }
  CHECK(kept_vessels_.empty());

  // Free old parts.  This must be done before binding the vessels, otherwise
  // the part subsets for the affected vessels will contain deleted parts.
  for (not_null<Vessel*> const vessel : loaded_vessels_) {
    vessel->FreeParts();
  }

  // Bind the vessels.  This guarantees that all part subsets are disjoint
  // unions of vessels.
  for (auto const& [_, vessel] : vessels_) {
    vessel->ForSomePart([&vessel = vessel](Part& first_part) {
      vessel->ForAllParts([&first_part](Part& part) {
        Subset<Part>::Unite(Subset<Part>::Find(first_part),
                            Subset<Part>::Find(part));
      });
    });
  }

  // Don't keep the grounded vessels.  This only destroys entire part subsets,
  // since being grounded is a subset property, and at this point part subsets
  // are disjoint unions of vessels.
  {
    // Note that we need to go through an intermediate set, since destroying a
    // vessel destroys its parts, which invalidates the intrusive `Subset` data
    // structure.
    VesselSet grounded_vessels;
    for (auto const& [_, vessel] : vessels_) {
      vessel->ForSomePart([&vessel = vessel, &grounded_vessels](Part& part) {
        if (Subset<Part>::Find(part).properties().grounded()) {
          grounded_vessels.insert(vessel.get());
        }
      });
    }
    for (not_null<Vessel*> const vessel : grounded_vessels) {
      loaded_vessels_.erase(vessel);
      LOG(INFO) << "Removing grounded vessel " << vessel->ShortDebugString();
      renderer_->ClearTargetVesselIf(vessel);
      zombie_prediction_adaptive_step_parameters_.insert_or_assign(
          vessel->guid(), vessel->prediction_adaptive_step_parameters());
      CHECK_EQ(vessels_.erase(vessel->guid()), 1);
    }
  }

  // We only need to collect one part per vessel, since the other parts are in
  // the same subset.
  for (auto const& [_, vessel] : vessels_) {
    Instant const vessel_time =
        is_loaded(vessel.get()) ? current_time_ - Δt : current_time_;
    vessel->ForSomePart([&vessel_time, this](Part& first_part) {
      Subset<Part>::Find(first_part).mutable_properties().Collect(
          pile_ups_,
          vessel_time,
          psychohistory_parameters_,
          history_fixed_step_parameters_,
          ephemeris_.get());
    });
  }

  // Now that the composition of the vessels is known, as well as their
  // intrinsic forces and torques, we may detect collapsibility changes.
  for (auto const& [_, vessel] : vessels_) {
    vessel->DetectCollapsibilityChange();
  }
}

void Plugin::SetPartApparentRigidMotion(
    PartId const part_id,
    RigidMotion<EccentricPart, ApparentWorld> const& rigid_motion) {
  // As a reference frame, `Apparent` differs from `World` only by having the
  // same axes as `Barycentric` and being nonrotating.  However, there is
  // another semantic distinction: `Apparent...` coordinates are uncorrected
  // data from the game, given immediately after its physics step; before using
  // them, we must correct them in accordance with the data computed by the pile
  // up.  This correction overrides the origin of position and velocity, so we
  // need not worry about the current definition of
  // `{World::origin, World::unmoving}` as we do when getting the actual degrees
  // of freedom (via `Plugin::BarycentricToWorld`).
  RigidMotion<ApparentWorld, Apparent> world_to_apparent{
      RigidTransformation<ApparentWorld, Apparent>{
          ApparentWorld::origin,
          Apparent::origin,
          OrthogonalMap<Barycentric, Apparent>::Identity() *
              renderer_->WorldToBarycentric(PlanetariumRotation()) *
              OrthogonalMap<ApparentWorld, World>::Identity()},
      Identity<Barycentric, Apparent>()(angular_velocity_of_world_),
      Apparent::unmoving};

  not_null<Vessel*> const vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  not_null<Part*> const part = vessel->part(part_id);
  CHECK(part->is_piled_up());
  part->containing_pile_up()->SetPartApparentRigidMotion(
      part,
      world_to_apparent * rigid_motion * part->MakeRigidToEccentricMotion());
}

RigidMotion<EccentricPart, World> Plugin::GetPartActualMotion(
    PartId const part_id,
    RigidMotion<Barycentric, World> const& barycentric_to_world) const {
  Part const& part = *FindOrDie(part_id_to_vessel_, part_id)->part(part_id);
  return barycentric_to_world * part.rigid_motion() *
         part.MakeRigidToEccentricMotion().Inverse();
}

DegreesOfFreedom<World> Plugin::CelestialWorldDegreesOfFreedom(
    Index const index,
    RigidMotion<Barycentric, World> const& barycentric_to_world,
    Instant const& time) const {
  return barycentric_to_world(
      FindOrDie(celestials_, index)->current_degrees_of_freedom(time));
}

RigidMotion<Barycentric, World> Plugin::BarycentricToWorld(
    bool const reference_part_is_unmoving,
    PartId const reference_part_id,
    std::optional<Position<World>> const& main_body_centre) const {
  BodyCentredNonRotatingReferenceFrame<Barycentric, MainBodyCentred> const
      main_body_frame{ephemeris_.get(), main_body_};

  auto const barycentric_to_main_body_motion =
      main_body_frame.ToThisFrameAtTime(current_time_);
  auto const barycentric_to_main_body_rotation =
      barycentric_to_main_body_motion.rigid_transformation().linear_map();
  Part const& reference_part = *FindOrDie(part_id_to_vessel_, reference_part_id)
                                    ->part(reference_part_id);
  auto const reference_part_degrees_of_freedom =
      barycentric_to_main_body_motion(reference_part.rigid_motion()(
          reference_part.MakeRigidToEccentricMotion().Inverse()(
              {EccentricPart::origin, EccentricPart::unmoving})));

  RigidTransformation<MainBodyCentred, World> const
      main_body_to_world_rigid_transformation = [&]() {
    if (main_body_centre.has_value()) {
      return RigidTransformation<World, MainBodyCentred>{
                 *main_body_centre,
                 MainBodyCentred::origin,
                 barycentric_to_main_body_rotation *
                     renderer_->WorldToBarycentric(
                         PlanetariumRotation())}.Inverse();
    } else {
      return RigidTransformation<MainBodyCentred, World>{
                 reference_part_degrees_of_freedom.position(),
                 World::origin,
                 renderer_->BarycentricToWorld(PlanetariumRotation()) *
                     barycentric_to_main_body_rotation.Inverse()};
    }
  }();
  RigidMotion<MainBodyCentred, World> const main_body_to_world = [&]() {
    if (reference_part_is_unmoving) {
      return RigidMotion<MainBodyCentred, World>{
          main_body_to_world_rigid_transformation,
          barycentric_to_main_body_rotation(angular_velocity_of_world_),
          reference_part_degrees_of_freedom.velocity()};
    } else {
      return RigidMotion<World, MainBodyCentred>{
          main_body_to_world_rigid_transformation.Inverse(),
          -(main_body_to_world_rigid_transformation.linear_map() *
                barycentric_to_main_body_rotation)(angular_velocity_of_world_),
          /*velocity_of_to_frame_origin=*/World::unmoving}.Inverse();
    }
  }();
  return main_body_to_world * barycentric_to_main_body_motion;
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  CHECK(!initializing_);
  CHECK_GT(t, current_time_);

  for (not_null<Vessel*> const vessel : loaded_vessels_) {
    vessel->ClearAllIntrinsicForcesAndTorques();
  }

  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
  ephemeris_->Prolong(current_time_).IgnoreError();
  UpdatePlanetariumRotation();
  loaded_vessels_.clear();
}

void Plugin::CatchUpLaggingVessels(VesselSet& collided_vessels) {
  CHECK(!initializing_);

  // Start all the integrations in parallel.
  std::vector<PileUpFuture> pile_up_futures;
  for (auto* const pile_up : pile_ups_) {
    pile_up_futures.emplace_back(
        pile_up,
        vessel_thread_pool_.Add([this, pile_up]() {
          // Note that there cannot be contention in the following method as
          // no two pile-ups are advanced at the same time.
          return pile_up->DeformAndAdvanceTime(current_time_);
        }));
  }

  // Wait for the integrations to finish and figure out which vessels collided
  // with a celestial.
  for (auto& pile_up_future : pile_up_futures) {
    WaitForVesselToCatchUp(pile_up_future, collided_vessels);
  }

  // Update the vessels.
  for (auto const& [_, vessel] : vessels_) {
    if (vessel->psychohistory()->back().time < current_time_) {
      if (Contains(collided_vessels, vessel.get())) {
        vessel->DisableDownsampling();
      }
      vessel->AdvanceTime();
    }
  }
}

not_null<std::unique_ptr<PileUpFuture>> Plugin::CatchUpVessel(
    GUID const& vessel_guid) {
  CHECK(!initializing_);

  // Find the vessel and the pile-up that contains it.
  Vessel& vessel = *FindOrDie(vessels_, vessel_guid);
  PileUp* pile_up = nullptr;
  vessel.ForSomePart([&pile_up](Part& part) {
    pile_up = part.containing_pile_up();
  });

  return make_not_null_unique<PileUpFuture>(
      pile_up,
      vessel_thread_pool_.Add([this, pile_up, &vessel]() {
        // Note that there can be contention in the following method if the
        // caller is catching-up two vessels belonging to the same pile-up in
        // parallel.
        absl::Status const status =
            pile_up->DeformAndAdvanceTime(current_time_);
        if (!status.ok()) {
          vessel.DisableDownsampling();
        }
        vessel.AdvanceTime();
        return status;
      }));
}

void Plugin::WaitForVesselToCatchUp(PileUpFuture& pile_up_future,
                                    VesselSet& collided_vessels) {
  PileUp const* const pile_up = pile_up_future.pile_up;
  auto& future = pile_up_future.future;
  future.wait();
  absl::Status const status = future.get();
  if (!status.ok()) {
    for (not_null<Part*> const part : pile_up->parts()) {
      not_null<Vessel*> const vessel =
          FindOrDie(part_id_to_vessel_, part->part_id());
      if (bool const inserted = collided_vessels.insert(vessel).second;
          inserted) {
        LOG(WARNING) << "Vessel " << vessel->ShortDebugString()
                     << " collided with a celestial: " << status.ToString();
      }
    }
  }
}

RelativeDegreesOfFreedom<AliceSun> Plugin::VesselFromParent(
    Index const parent_index,
    GUID const& vessel_guid) const {
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      FindOrDie(vessels_, vessel_guid);
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  if (vessel->parent() != parent) {
    vessel->set_parent(parent);
  }
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      vessel->psychohistory()->back().degrees_of_freedom -
      vessel->parent()->current_degrees_of_freedom(current_time_);
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  return result;
}

RelativeDegreesOfFreedom<AliceSun> Plugin::CelestialFromParent(
    Index const celestial_index) const {
  CHECK(!initializing_);
  ephemeris_->Prolong(current_time_).IgnoreError();
  Celestial const& celestial = *FindOrDie(celestials_, celestial_index);
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      celestial.current_degrees_of_freedom(current_time_) -
      celestial.parent()->current_degrees_of_freedom(current_time_);
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  return result;
}

void Plugin::SetPredictionAdaptiveStepParameters(
    GUID const& vessel_guid,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        prediction_adaptive_step_parameters) const {
  // If there is a target vessel, it is integrated with the same parameters as
  // the given (current) vessel.  This makes it possible to plot the prediction
  // of the current vessel.
  auto& vessel = *FindOrDie(vessels_, vessel_guid);
  if (renderer_->HasTargetVessel()) {
    auto& target_vessel = renderer_->GetTargetVessel();
    if (vessel.has_flight_plan()) {
      // If there is a target vessel and our vessel has a flight plan, make sure
      // that we don't shorten the target vessel's prediction as it would hurt
      // our flight plan.
      auto const target_prediction_adaptive_step_parameters =
        target_vessel.prediction_adaptive_step_parameters();
      target_vessel.set_prediction_adaptive_step_parameters(
          Ephemeris<Barycentric>::AdaptiveStepParameters(
              prediction_adaptive_step_parameters.integrator(),
              std::max(target_prediction_adaptive_step_parameters.max_steps(),
                       prediction_adaptive_step_parameters.max_steps()),
              prediction_adaptive_step_parameters
                  .length_integration_tolerance(),
              prediction_adaptive_step_parameters
                  .speed_integration_tolerance()));
    } else {
      target_vessel.set_prediction_adaptive_step_parameters(
          prediction_adaptive_step_parameters);
    }
  }
  vessel.set_prediction_adaptive_step_parameters(
      prediction_adaptive_step_parameters);
}

void Plugin::UpdatePrediction(std::vector<GUID> const& vessel_guids) const {
  CHECK(!initializing_);
  std::set<not_null<Vessel*>> predicted_vessels;
  for (auto const& guid : vessel_guids) {
    predicted_vessels.insert(FindOrDie(vessels_, guid).get());
  }
  Vessel* target_vessel = nullptr;

  // If there is a target vessel, ensure that the prediction of the
  // `predicted_vessels` is not longer than that of the target vessel.  This is
  // necessary to build the targeting frame.
  if (renderer_->HasTargetVessel()) {
    target_vessel = &renderer_->GetTargetVessel();
    target_vessel->RefreshPrediction();
    for (auto const vessel : predicted_vessels) {
      vessel->RefreshPrediction(target_vessel->prediction()->back().time);
    }
  } else {
    for (auto const vessel : predicted_vessels) {
      vessel->RefreshPrediction();
    }
  }
  for (auto const& [guid, vessel] : vessels_) {
    if (!Contains(predicted_vessels, vessel.get()) &&
        vessel.get() != target_vessel) {
      vessel->StopPrognosticator();
    }
  }
}

void Plugin::CreateFlightPlan(GUID const& vessel_guid,
                              Instant const& final_time,
                              Mass const& initial_mass) const {
  CHECK(!initializing_);
  // TODO(phl): Serialize the burn parameters.  We should also probably
  // distinguish the coast parameters from the prediction parameters.
  FindOrDie(vessels_, vessel_guid)->CreateFlightPlan(
      final_time,
      initial_mass,
      DefaultPredictionParameters(),
      DefaultBurnParameters());
}

void Plugin::ExtendPredictionForFlightPlan(GUID const& vessel_guid) const {
  auto& vessel = *FindOrDie(vessels_, vessel_guid);

  // If there is a target vessel, and the prediction of the target is too short
  // for the flight plan, multiply the number of steps by 4 (one notch in the
  // UI).  We don't wait for the prediction to be computed, though, so there is
  // no guarantee that it will be long enough.  Presumably the user will keep
  // increasing the flight plan length and get what they want, ultimately.
  if (renderer_->HasTargetVessel() && vessel.has_flight_plan()) {
    vessel.ReadFlightPlanFromMessage();
    auto& target_vessel = renderer_->GetTargetVessel();
    if (target_vessel.prediction()->back().time <
        vessel.flight_plan().actual_final_time()) {
      auto prediction_adaptive_step_parameters =
          target_vessel.prediction_adaptive_step_parameters();
      prediction_adaptive_step_parameters.set_max_steps(
          std::min(max_steps_in_prediction,
                   prediction_adaptive_step_parameters.max_steps() * 4));
      target_vessel.set_prediction_adaptive_step_parameters(
          prediction_adaptive_step_parameters);
    }
  }
}

void Plugin::ComputeAndRenderApsides(
    Index const celestial_index,
    Trajectory<Barycentric> const& trajectory,
    DiscreteTrajectory<Barycentric>::iterator const& begin,
    DiscreteTrajectory<Barycentric>::iterator const& end,
    Instant const& t_max,
    Position<World> const& sun_world_position,
    int const max_points,
    DiscreteTrajectory<World>& apoapsides,
    DiscreteTrajectory<World>& periapsides) const {
  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(FindOrDie(celestials_, celestial_index)->trajectory(),
                 trajectory,
                 begin, end,
                 t_max,
                 max_points,
                 apoapsides_trajectory,
                 periapsides_trajectory);
  apoapsides = renderer_->RenderBarycentricTrajectoryInWorld(
                   current_time_,
                   apoapsides_trajectory.begin(),
                   apoapsides_trajectory.end(),
                   sun_world_position,
                   PlanetariumRotation());
  periapsides = renderer_->RenderBarycentricTrajectoryInWorld(
                    current_time_,
                    periapsides_trajectory.begin(),
                    periapsides_trajectory.end(),
                    sun_world_position,
                    PlanetariumRotation());
}

std::optional<DiscreteTrajectory<World>::value_type>
Plugin::ComputeAndRenderFirstCollision(
    Index const celestial_index,
    Trajectory<Barycentric> const& trajectory,
    DiscreteTrajectory<Barycentric>::iterator const& begin,
    DiscreteTrajectory<Barycentric>::iterator const& end,
    Position<World> const& sun_world_position,
    int max_points,
    std::function<Length(Angle const& latitude,
                         Angle const& longitude)> const& radius) const {
  auto const& celestial = FindOrDie(celestials_, celestial_index);
  auto const& celestial_body = *celestial->body();
  auto const& celestial_trajectory = celestial->trajectory();

  // TODO(phl): We should cache the apsides.
  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(celestial_trajectory,
                 trajectory,
                 begin, end,
                 /*t_max=*/InfiniteFuture,
                 max_points,
                 apoapsides_trajectory,
                 periapsides_trajectory);

  const auto intervals = ComputeCollisionIntervals(celestial_body,
                                                   celestial_trajectory,
                                                   trajectory,
                                                   apoapsides_trajectory,
                                                   periapsides_trajectory);

  VLOG(1) << "Found " << intervals.size() << " collision intervals";
  for (auto const& interval : intervals) {
    VLOG(1) << "Collision interval: " << interval;
    auto const maybe_collision = ComputeFirstCollision(celestial_body,
                                                       celestial_trajectory,
                                                       trajectory,
                                                       interval,
                                                       MaxCollisionError(),
                                                       radius);
    if (maybe_collision.has_value()) {
      auto const& collision = maybe_collision.value();

      // We create a trajectory with a single point to simplify rendering.
      DiscreteTrajectory<Barycentric> trajectory_to_render;
      CHECK_OK(trajectory_to_render.Append(collision.time,
                                           collision.degrees_of_freedom));
      DiscreteTrajectory<World> rendered_trajectory =
          renderer_->RenderBarycentricTrajectoryInWorld(
              current_time_,
              trajectory_to_render.begin(),
              trajectory_to_render.end(),
              sun_world_position,
              PlanetariumRotation());
      return rendered_trajectory.front();
    }
  }

  // No collision.
  return std::nullopt;
}

void Plugin::ComputeAndRenderClosestApproaches(
    Trajectory<Barycentric> const& trajectory,
    DiscreteTrajectory<Barycentric>::iterator const& begin,
    DiscreteTrajectory<Barycentric>::iterator const& end,
    Position<World> const& sun_world_position,
    int const max_points,
    DiscreteTrajectory<World>& closest_approaches) const {
  CHECK(renderer_->HasTargetVessel());

  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(*renderer_->GetTargetVessel().prediction(),
                 trajectory,
                 begin, end,
                 /*t_max=*/InfiniteFuture,
                 max_points,
                 apoapsides_trajectory,
                 periapsides_trajectory);
  closest_approaches =
      renderer_->RenderBarycentricTrajectoryInWorld(
          current_time_,
          periapsides_trajectory.begin(),
          periapsides_trajectory.end(),
          sun_world_position,
          PlanetariumRotation());
}

void Plugin::ComputeAndRenderNodes(
    DiscreteTrajectory<Barycentric>::iterator const& begin,
    DiscreteTrajectory<Barycentric>::iterator const& end,
    Instant const& t_max,
    Position<World> const& sun_world_position,
    int const max_points,
    std::vector<Renderer::Node>& ascending,
    std::vector<Renderer::Node>& descending) const {
  auto const trajectory_in_plotting =
      renderer_->RenderBarycentricTrajectoryInPlotting(begin, end);

  auto const* const cast_plotting_frame = dynamic_cast<
      BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation> const*>(
      &*renderer_->GetPlottingFrame());
  // The body-centred non rotating frame does not rotate; its reference plane is
  // not an inherent dynamical property.  When using this frame, discard the
  // nodes if they are far enough from the central body that its equator is
  // irrelevant.
  Length const threshold =
      cast_plotting_frame == nullptr
          ? Infinity<Length>
          : EquatorRelevanceThreshold(
                *dynamic_cast_not_null<RotatingBody<Barycentric> const*>(
                    cast_plotting_frame->centre()));
  auto const show_node = [threshold](DegreesOfFreedom<Navigation> const& dof) {
    return (dof.position() - Navigation::origin).Norm() < threshold;
  };

  DiscreteTrajectory<Navigation> ascending_trajectory;
  DiscreteTrajectory<Navigation> descending_trajectory;
  // The so-called North is orthogonal to the plane of the trajectory.
  ComputeNodes(trajectory_in_plotting,
               trajectory_in_plotting.begin(),
               trajectory_in_plotting.end(),
               t_max,
               Vector<double, Navigation>({0, 0, 1}),
               max_points,
               ascending_trajectory,
               descending_trajectory,
               show_node).IgnoreError();

  ascending = renderer_->RenderNodes(current_time_,
                                     ascending_trajectory.begin(),
                                     ascending_trajectory.end(),
                                     sun_world_position,
                                     PlanetariumRotation());
  descending = renderer_->RenderNodes(current_time_,
                                      descending_trajectory.begin(),
                                      descending_trajectory.end(),
                                      sun_world_position,
                                      PlanetariumRotation());
}

bool Plugin::HasCelestial(Index const index) const {
  return Contains(celestials_, index);
}

Celestial const& Plugin::GetCelestial(Index const index) const {
  return *FindOrDie(celestials_, index);
}

bool Plugin::HasVessel(GUID const& vessel_guid) const {
  return Contains(vessels_, vessel_guid);
}

not_null<Vessel*> Plugin::GetVessel(GUID const& vessel_guid) const {
  CHECK(!initializing_);
  return FindOrDie(vessels_, vessel_guid).get();
}

void Plugin::ClearOrbitAnalysersOfVesselsOtherThan(Vessel const& vessel) {
  for (auto const& [guid, v] : vessels_) {
    if (v.get() != &vessel) {
      v->ClearOrbitAnalyser();
    }
  }
}

not_null<std::unique_ptr<Planetarium>> Plugin::NewPlanetarium(
    Planetarium::Parameters const& parameters,
    Perspective<Navigation, Camera> const& perspective,
    Planetarium::PlottingToScaledSpaceConversion plotting_to_scaled_space)
    const {
  return make_not_null_unique<Planetarium>(parameters,
                                           perspective,
                                           ephemeris_.get(),
                                           renderer_->GetPlottingFrame(),
                                           std::move(plotting_to_scaled_space));
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBarycentricRotatingNavigationFrame(
    Index const primary_index,
    Index const secondary_index) const {
  CHECK(!initializing_);
  Celestial const& primary = *FindOrDie(celestials_, primary_index);
  Celestial const& secondary = *FindOrDie(celestials_, secondary_index);
  return make_not_null_unique<
      BarycentricRotatingReferenceFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          primary.body(),
          secondary.body());
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBodyCentredBodyDirectionNavigationFrame(
    Index const primary_index,
    Index const secondary_index) const {
  CHECK(!initializing_);
  Celestial const& primary = *FindOrDie(celestials_, primary_index);
  Celestial const& secondary = *FindOrDie(celestials_, secondary_index);
  return make_not_null_unique<
      BodyCentredBodyDirectionReferenceFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          primary.body(),
          secondary.body());
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBodyCentredNonRotatingNavigationFrame(
    Index const reference_body_index) const {
  CHECK(!initializing_);
  Celestial const& reference_body =
      *FindOrDie(celestials_, reference_body_index);
  return make_not_null_unique<
      BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          reference_body.body());
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBodySurfaceNavigationFrame(
    Index const reference_body_index) const {
  CHECK(!initializing_);
  Celestial const& reference_body =
      *FindOrDie(celestials_, reference_body_index);
  return make_not_null_unique<
      BodySurfaceReferenceFrame<Barycentric, Navigation>>(
      ephemeris_.get(), reference_body.body());
}

not_null<std::unique_ptr<PlottingFrame>>
Plugin::NewRotatingPulsatingPlottingFrame(
    std::vector<Index> const& primary_indices,
    std::vector<Index> const& secondary_indices) const {
  std::vector<not_null<MassiveBody const*>> primaries;
  for (Index const i : primary_indices) {
    Celestial const& primary = *FindOrDie(celestials_, i);
    primaries.push_back(primary.body());
  }
  std::vector<not_null<MassiveBody const*>> secondaries;
  for (Index const i : secondary_indices) {
    Celestial const& secondary = *FindOrDie(celestials_, i);
    secondaries.push_back(secondary.body());
  }
  return make_not_null_unique<
      RotatingPulsatingReferenceFrame<Barycentric, Navigation>>(
      ephemeris_.get(), std::move(primaries), std::move(secondaries));
}

void Plugin::SetTargetVessel(GUID const& vessel_guid,
                             Index const reference_body_index) {
  not_null<Celestial const*> const celestial =
      FindOrDie(celestials_, reference_body_index).get();
  not_null<Vessel*> const vessel = FindOrDie(vessels_, vessel_guid).get();
  renderer_->SetTargetVessel(vessel, celestial, ephemeris_.get());
}

std::unique_ptr<FrameField<World, Navball>> Plugin::NavballFrameField(
    Position<World> const& sun_world_position) const {

  using RightHandedNavball = Frame<struct RightHandedNavballTag>;

  // TODO(phl): Clean up this mess!
  class NavballFrameField : public FrameField<World, Navball> {
   public:
    NavballFrameField(
        not_null<Plugin const*> const plugin,
        not_null<std::unique_ptr<FrameField<Navigation, RightHandedNavball>>>
            navigation_right_handed_field,
        Position<World> const& sun_world_position)
        : plugin_(plugin),
          navigation_right_handed_field_(
              std::move(navigation_right_handed_field)),
          sun_world_position_(sun_world_position) {}

    NavballFrameField(
        not_null<Plugin const*> const plugin,
        not_null<std::unique_ptr<FrameField<Barycentric, RightHandedNavball>>>
            barycentric_right_handed_field,
        Position<World> const& sun_world_position)
        : plugin_(plugin),
          barycentric_right_handed_field_(
              std::move(barycentric_right_handed_field)),
          sun_world_position_(sun_world_position) {}

    Rotation<Navball, World> FromThisFrame(
        Position<World> const& q) const override {
      Instant const& current_time = plugin_->current_time_;
      auto const& planetarium_rotation = plugin_->PlanetariumRotation();
      auto const& renderer = plugin_->renderer();
      plugin_->ephemeris_->Prolong(current_time).IgnoreError();

      Position<Navigation> const q_in_plotting =
          renderer.WorldToPlotting(current_time,
                                   sun_world_position_,
                                   planetarium_rotation)(q);

      OrthogonalMap<RightHandedNavball,
                    Barycentric> const right_handed_navball_to_barycentric =
          barycentric_right_handed_field_ == nullptr
              ? renderer.PlottingToBarycentric(current_time)
                        .orthogonal_map¹₁() *
                    navigation_right_handed_field_->FromThisFrame(q_in_plotting)
                        .Forget<OrthogonalMap>()
              : barycentric_right_handed_field_
                    ->FromThisFrame(
                        renderer.WorldToBarycentric(current_time,
                                                    sun_world_position_,
                                                    planetarium_rotation)(q))
                    .Forget<OrthogonalMap>();

      // KSP's navball has x west, y up, z south.
      // We want x north, y east, z down.
      OrthogonalMap<Navball, World> const orthogonal_map =
          renderer.BarycentricToWorld(planetarium_rotation) *
          right_handed_navball_to_barycentric *
          Permutation<World, RightHandedNavball>(OddPermutation::XZY)
              .Forget<OrthogonalMap>() *
          Rotation<Navball, World>(π / 2 * Radian,
                                   Bivector<double, World>({0, 1, 0}),
                                   DefinesFrame<Navball>())
                                   .Forget<OrthogonalMap>();
      return orthogonal_map.AsRotation();
    }

   private:
    not_null<Plugin const*> const plugin_;
    std::unique_ptr<FrameField<Navigation, RightHandedNavball>> const
        navigation_right_handed_field_;
    std::unique_ptr<FrameField<Barycentric, RightHandedNavball>> const
        barycentric_right_handed_field_;
    Position<World> const sun_world_position_;
  };

  std::unique_ptr<FrameField<Navigation, RightHandedNavball>> frame_field;
  auto* const plotting_frame_as_body_surface_reference_frame =
      dynamic_cast<BodySurfaceReferenceFrame<Barycentric, Navigation> const*>(
          &*renderer_->GetPlottingFrame());
  if (plotting_frame_as_body_surface_reference_frame == nullptr) {
    return std::make_unique<NavballFrameField>(
        this,
        make_not_null_unique<
            CoordinateFrameField<Navigation, RightHandedNavball>>(),
        sun_world_position);
  } else {
    return std::make_unique<NavballFrameField>(
        this,
        make_not_null_unique<
            BodySurfaceFrameField<Barycentric, RightHandedNavball>>(
                *ephemeris_,
                current_time_,
                plotting_frame_as_body_surface_reference_frame->centre()),
        sun_world_position);
  }
}

Vector<double, World> Plugin::VesselTangent(GUID const& vessel_guid) const {
  return renderer_->FrenetToWorld(
      *FindOrDie(vessels_, vessel_guid),
      PlanetariumRotation())(Vector<double, Frenet<Navigation>>({1, 0, 0}));
}

Vector<double, World> Plugin::VesselNormal(GUID const& vessel_guid) const {
  return renderer_->FrenetToWorld(
      *FindOrDie(vessels_, vessel_guid),
      PlanetariumRotation())(Vector<double, Frenet<Navigation>>({0, 1, 0}));
}

Vector<double, World> Plugin::VesselBinormal(GUID const& vessel_guid) const {
  return renderer_->FrenetToWorld(
      *FindOrDie(vessels_, vessel_guid),
      PlanetariumRotation())(Vector<double, Frenet<Navigation>>({0, 0, 1}));
}

Velocity<World> Plugin::UnmanageableVesselVelocity(
    RelativeDegreesOfFreedom<AliceSun> const& degrees_of_freedom,
    Index const parent_index) const {
  auto const parent_degrees_of_freedom =
      FindOrDie(celestials_,
                parent_index)->current_degrees_of_freedom(current_time_);
  return VesselVelocity(
      current_time_,
      parent_degrees_of_freedom +
          PlanetariumRotation().Inverse()(degrees_of_freedom));
}

Velocity<World> Plugin::VesselVelocity(GUID const& vessel_guid) const {
  Vessel const& vessel = *FindOrDie(vessels_, vessel_guid);
  auto const back = vessel.psychohistory()->back();
  return VesselVelocity(back.time, back.degrees_of_freedom);
}

void Plugin::RequestReanimation(Instant const& desired_t_min) const {
  ephemeris_->RequestReanimation(desired_t_min);
}

Instant Plugin::GameEpoch() const {
  return game_epoch_;
}

Instant Plugin::CurrentTime() const {
  return current_time_;
}

Rotation<Barycentric, AliceSun> const& Plugin::PlanetariumRotation() const {
  return *cached_planetarium_rotation_;
}

Rotation<CameraCompensatedReference, CameraReference> const&
Plugin::CameraCompensation() const {
  return *camera_compensation_;
}

Renderer& Plugin::renderer() {
  return *renderer_;
}

Renderer const& Plugin::renderer() const {
  return *renderer_;
}

GeometricPotentialPlotter& Plugin::geometric_potential_plotter() {
  CHECK(!initializing_);
  return *geometric_potential_plotter_;
}

GeometricPotentialPlotter const& Plugin::geometric_potential_plotter() const {
  CHECK(!initializing_);
  return *geometric_potential_plotter_;
}

void Plugin::WriteToMessage(
    not_null<serialization::Plugin*> const message) const {
  LOG(INFO) << __FUNCTION__;
  CHECK(!initializing_);
  message->set_uses_correct_sin_cos(uses_correct_sin_cos_);
  if (system_fingerprint_ != 0) {
    message->set_system_fingerprint(system_fingerprint_);
  }
  ephemeris_->Prolong(current_time_).IgnoreError();
  std::map<not_null<Celestial const*>, Index const> celestial_to_index;
  for (auto const& [index, owned_celestial] : celestials_) {
    celestial_to_index.emplace(owned_celestial.get(), index);
  }
  for (auto const& [index, owned_celestial] : celestials_) {
    auto* const celestial_message = message->add_celestial();
    celestial_message->set_index(index);
    if (owned_celestial->has_parent()) {
      Index const parent_index =
          FindOrDie(celestial_to_index, owned_celestial->parent());
      celestial_message->set_parent_index(parent_index);
    }
    celestial_message->set_ephemeris_index(
        ephemeris_->serialization_index_for_body(owned_celestial->body()));
  }

  // Construct a map to help serialization of the pile-ups.
  std::map<not_null<PileUp const*>, int> serialization_index_to_pile_up;
  int serialization_index = 0;
  for (auto const* pile_up : pile_ups_) {
    serialization_index_to_pile_up[pile_up] = serialization_index++;
  }
  auto const serialization_index_for_pile_up =
      [&serialization_index_to_pile_up](not_null<PileUp const*> const pile_up) {
        return serialization_index_to_pile_up.at(pile_up);
      };

  std::map<not_null<Vessel const*>, GUID const> vessel_to_guid;
  for (auto const& [guid, vessel] : vessels_) {
    vessel_to_guid.emplace(vessel.get(), guid);
    auto* const vessel_message = message->add_vessel();
    vessel_message->set_guid(guid);
    vessel->WriteToMessage(vessel_message->mutable_vessel(),
                           serialization_index_for_pile_up);
    Index const parent_index = FindOrDie(celestial_to_index, vessel->parent());
    vessel_message->set_parent_index(parent_index);
    vessel_message->set_loaded(Contains(loaded_vessels_, vessel.get()));
    vessel_message->set_kept(Contains(kept_vessels_, vessel.get()));
  }
  for (auto const& [part_id, vessel] : part_id_to_vessel_) {
    (*message->mutable_part_id_to_vessel())[part_id] = vessel_to_guid[vessel];
  }
  for (auto const& [guid, parameters] :
       zombie_prediction_adaptive_step_parameters_) {
    auto* const zombie_message = message->add_zombie();
    zombie_message->set_guid(guid);
    parameters.WriteToMessage(zombie_message->mutable_prediction_parameters());
  }

  ephemeris_->WriteToMessage(message->mutable_ephemeris());

  // `history_downsampling_parameters_` is not persisted.
  history_fixed_step_parameters_.WriteToMessage(
      message->mutable_history_parameters());
  psychohistory_parameters_.WriteToMessage(
      message->mutable_psychohistory_parameters());

  planetarium_rotation_.WriteToMessage(message->mutable_planetarium_rotation());
  game_epoch_.WriteToMessage(message->mutable_game_epoch());
  current_time_.WriteToMessage(message->mutable_current_time());
  Index const sun_index = FindOrDie(celestial_to_index, sun_);
  message->set_sun_index(sun_index);
  renderer_->WriteToMessage(message->mutable_renderer());

  for (auto* const pile_up : pile_ups_) {
    pile_up->WriteToMessage(message->add_pile_up());
  }
}

not_null<std::unique_ptr<Plugin>> Plugin::ReadFromMessage(
    serialization::Plugin const& message) {
  LOG(INFO) << __FUNCTION__;

  auto const history_parameters =
      Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
          message.history_parameters());
  auto const psychohistory_parameters =
      Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
          message.psychohistory_parameters());
  not_null<std::unique_ptr<Plugin>> plugin =
      std::unique_ptr<Plugin>(new Plugin(history_parameters,
                                         psychohistory_parameters));

  plugin->uses_correct_sin_cos_ = message.has_uses_correct_sin_cos() &&
      message.uses_correct_sin_cos();

  if (message.has_system_fingerprint()) {
    plugin->system_fingerprint_ = message.system_fingerprint();
    std::string details = "this is an unknown system";
    for (auto const ksp_version : {KSP122, KSP191}) {
      if (plugin->system_fingerprint_ ==
          KSPStockSystemFingerprints[ksp_version]) {
        details = "this is the dreaded KSP stock system!";
        break;
      } else if (plugin->system_fingerprint_ ==
                 KSPStabilizedSystemFingerprints[ksp_version]) {
        details = "this is the stabilized KSP system, all hail retrobop!";
        break;
      }
    }
    LOG(INFO) << "System has fingerprint 0x" << std::hex << std::uppercase
              << plugin->system_fingerprint_ << "; " << details;
  }

  plugin->game_epoch_ = Instant::ReadFromMessage(message.game_epoch());
  plugin->current_time_ = Instant::ReadFromMessage(message.current_time());
  plugin->planetarium_rotation_ =
      Angle::ReadFromMessage(message.planetarium_rotation());

  // The ephemeris constructed here is *not* prolonged and needs to be
  // explicitly prolonged to cover all the instants that we care about.  Note
  // that it is important to pick the most recent checkpoint that covers the
  // current time: an older checkpoint would require unnecessary work in
  // Prolong that could be postponed until reanimation; a newer checkpoint would
  // not cover the current time.
  plugin->ephemeris_ =
      Ephemeris<Barycentric>::ReadFromMessage(/*using_checkpoint_at_or_before=*/
                                              plugin->current_time_,
                                              message.ephemeris());
  plugin->ephemeris_->Prolong(plugin->game_epoch_).IgnoreError();
  plugin->ephemeris_->Prolong(plugin->current_time_).IgnoreError();
  CHECK_LE(plugin->ephemeris_->t_min(), plugin->current_time_);

  ReadCelestialsFromMessages(*plugin->ephemeris_,
                             message.celestial(),
                             plugin->celestials_,
                             plugin->name_to_index_);

  for (auto const& vessel_message : message.vessel()) {
    not_null<Celestial const*> const parent =
        FindOrDie(plugin->celestials_, vessel_message.parent_index()).get();
    not_null<std::unique_ptr<Vessel>> vessel = Vessel::ReadFromMessage(
        vessel_message.vessel(),
        parent,
        plugin->ephemeris_.get(),
        [&part_id_to_vessel = plugin->part_id_to_vessel_](
            PartId const part_id) {
          CHECK_NE(part_id_to_vessel.erase(part_id), 0) << part_id;
        });

    if (vessel_message.loaded()) {
      plugin->loaded_vessels_.insert(vessel.get());
    }
    if (vessel_message.kept()) {
      plugin->kept_vessels_.insert(vessel.get());
    }
    bool const inserted = plugin->vessels_.emplace(
        vessel_message.guid(), std::move(vessel)).second;
    CHECK(inserted);
  }

  for (auto const& [part_id, guid] : message.part_id_to_vessel()) {
    auto const& vessel = FindOrDie(plugin->vessels_, guid);
    plugin->part_id_to_vessel_.emplace(part_id, vessel.get());
  }

  plugin->sun_ = FindOrDie(plugin->celestials_, message.sun_index()).get();
  plugin->main_body_ = plugin->sun_->body();
  plugin->UpdatePlanetariumRotation();

  bool const is_pre_cauchy = message.has_pre_cauchy_plotting_frame();
  LOG_IF(WARNING, is_pre_cauchy) << "Reading pre-Cauchy Plugin";

  if (is_pre_cauchy) {
    plugin->renderer_ =
        std::make_unique<Renderer>(
            plugin->sun_,
            NavigationFrame::ReadFromMessage(
                message.pre_cauchy_plotting_frame(),
                plugin->ephemeris_.get()));
  } else {
    plugin->renderer_ = Renderer::ReadFromMessage(message.renderer(),
                                                  plugin->sun_,
                                                  plugin->ephemeris_.get());
  }

  plugin->geometric_potential_plotter_.emplace(plugin->ephemeris_.get());

  // Note that for proper deserialization of parts this list must be
  // reconstructed in its original order.
  auto const part_id_to_part =
      [&part_id_to_vessel =
            plugin->part_id_to_vessel_](PartId const part_id) {
        not_null<Vessel*> const vessel = part_id_to_vessel.at(part_id);
        not_null<Part*> const part = vessel->part(part_id);
        return part;
      };
  for (auto const& pile_up_message : message.pile_up()) {
    // First push a nullptr to be able to capture an iterator to the new
    // location in the list in the deletion callback.
    plugin->pile_ups_.push_back(nullptr);
    auto deletion_callback = [it = std::prev(plugin->pile_ups_.end()),
                              &pile_ups = plugin->pile_ups_]() {
      pile_ups.erase(it);
    };
    auto const pile_up = PileUp::ReadFromMessage(pile_up_message,
                                                 part_id_to_part,
                                                 plugin->ephemeris_.get(),
                                                 std::move(deletion_callback))
                             .release();
    *plugin->pile_ups_.rbegin() = pile_up;
  }

  // Now fill the containing pile-up of all the parts.  This gives ownership of
  // the pile-ups to the parts.  To do that, we first build shared pointers for
  // all the pile-ups.
  std::vector<not_null<std::shared_ptr<PileUp>>> shared_pile_ups;
  for (auto* const pile_up : plugin->pile_ups_) {
    shared_pile_ups.emplace_back(check_not_null(pile_up));
  }
  auto const pile_up_for_serialization_index =
      [&shared_pile_ups](int const serialization_index) {
    return shared_pile_ups.at(serialization_index);
  };

  for (auto const& vessel_message : message.vessel()) {
    GUID const guid = vessel_message.guid();
    auto const& vessel = FindOrDie(plugin->vessels_, guid);
    vessel->FillContainingPileUpsFromMessage(vessel_message.vessel(),
                                             pile_up_for_serialization_index);
  }

  plugin->initializing_.Flop();
  return plugin;
}

Plugin::Plugin(
    Ephemeris<Barycentric>::FixedStepParameters history_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters
        psychohistory_parameters)
    : history_downsampling_parameters_(DefaultDownsamplingParameters()),
      history_fixed_step_parameters_(std::move(history_parameters)),
      psychohistory_parameters_(std::move(psychohistory_parameters)),
      vessel_thread_pool_(
          /*pool_size=*/2 * std::thread::hardware_concurrency()) {}

void Plugin::InitializeIndices(std::string const& name,
                               Index const celestial_index,
                               std::optional<Index> const& parent_index) {
  bool inserted = name_to_index_.emplace(name, celestial_index).second;
  CHECK(inserted) << name;
  inserted = index_to_name_.emplace(celestial_index, name).second;
  CHECK(inserted) << celestial_index;
  inserted = parents_.emplace(celestial_index, parent_index).second;
  CHECK(inserted) << celestial_index;
}

void Plugin::UpdatePlanetariumRotation() {
  using PlanetariumFrame = Frame<struct PlanetariumFrameTag>;

  CHECK_NOTNULL(main_body_);
  Rotation<Barycentric, PlanetariumFrame> const to_planetarium =
      main_body_->ToCelestialFrame<PlanetariumFrame>();
  cached_planetarium_rotation_ =
      Rotation<PlanetariumFrame, AliceSun>(
          planetarium_rotation_,
          Bivector<double, PlanetariumFrame>({0, 0, 1}),
          DefinesFrame<AliceSun>{}) *
      to_planetarium;
  camera_compensation_ = Rotation<CameraCompensatedReference, CameraReference>(
      -planetarium_rotation_,
      Bivector<double, CameraReference>({0, 1, 0}),
      DefinesFrame<CameraCompensatedReference>{});
}

Velocity<World> Plugin::VesselVelocity(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) const {
  DegreesOfFreedom<Navigation> const plotting_frame_degrees_of_freedom =
      renderer_->BarycentricToPlotting(time)(degrees_of_freedom);
  // Note that in the rotating-pulsating reference frame, this value is given in
  // current metres per second (that is, in metres at `time` per second, for a
  // velocity at `time`).
  return renderer_->PlottingToWorld(time, PlanetariumRotation())(
      plotting_frame_degrees_of_freedom.velocity());
}

template<typename T>
void Plugin::ReadCelestialsFromMessages(
    Ephemeris<Barycentric> const& ephemeris,
    google::protobuf::RepeatedPtrField<T> const& celestial_messages,
    IndexToOwnedCelestial& celestials,
    std::map<std::string, Index>& name_to_index) {
  auto const& bodies = ephemeris.bodies();
  int index = 0;
  for (auto const& celestial_message : celestial_messages) {
    bool const is_pre_cauchy = !celestial_message.has_ephemeris_index();
    LOG_IF(WARNING, is_pre_cauchy) << "Reading pre-Cauchy Celestial";

    auto const& body = is_pre_cauchy
                           ? bodies[index++]
                           : bodies[celestial_message.ephemeris_index()];
    auto [it, inserted] = celestials.emplace(
        celestial_message.index(),
        make_not_null_unique<Celestial>(
            dynamic_cast_not_null<RotatingBody<Barycentric> const*>(
                body)));
    CHECK(inserted) << celestial_message.index();
    it->second->set_trajectory(ephemeris.trajectory(body));

    inserted =
        name_to_index.emplace(body->name(), celestial_message.index()).second;
    CHECK(inserted) << body->name();
  }
  for (auto const& celestial_message : celestial_messages) {
    if (celestial_message.has_parent_index()) {
      not_null<std::unique_ptr<Celestial>> const& celestial =
          FindOrDie(celestials, celestial_message.index());
      not_null<Celestial const*> const parent =
          FindOrDie(celestials, celestial_message.parent_index()).get();
      celestial->set_parent(parent);
    }
  }
}

template<typename... Args>
void Plugin::AddPart(not_null<Vessel*> const vessel,
                     PartId const part_id,
                     std::string const& name,
                     Args... args) {
  auto const [_, inserted] = part_id_to_vessel_.emplace(part_id, vessel);
  CHECK(inserted) << NAMED(part_id);
  auto deletion_callback = [part_id, &map = part_id_to_vessel_] {
    // This entails a lookup, but iterators are not stable in `flat_hash_map`.
    map.erase(part_id);
  };
  auto part = make_not_null_unique<Part>(part_id,
                                         name,
                                         std::forward<Args>(args)...,
                                         std::move(deletion_callback));
  vessel->AddPart(std::move(part));
}

bool Plugin::is_loaded(not_null<Vessel*> vessel) const {
  return Contains(loaded_vessels_, vessel);
}

}  // namespace internal
}  // namespace _plugin
}  // namespace ksp_plugin
}  // namespace principia
