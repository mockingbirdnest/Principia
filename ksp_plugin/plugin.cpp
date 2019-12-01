
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <list>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/solar_system_fingerprints.hpp"
#include "astronomy/stabilize_ksp.hpp"
#include "astronomy/time_scales.hpp"
#include "base/file.hpp"
#include "base/hexadecimal.hpp"
#include "base/map_util.hpp"
#include "base/not_null.hpp"
#include "base/optional_logging.hpp"
#include "base/serialization.hpp"
#include "base/status.hpp"
#include "base/unique_ptr_logging.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/identity.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/permutation.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "glog/logging.h"
#include "glog/stl_logging.h"
#include "ksp_plugin/equator_relevance_threshold.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/part_subsets.hpp"
#include "physics/apsides.hpp"
#include "physics/barycentric_rotating_dynamic_frame_body.hpp"
#include "physics/body_centred_body_direction_dynamic_frame.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/body_surface_frame_field.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/frame_field.hpp"
#include "physics/massive_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using astronomy::KSPStabilizedSystemFingerprint;
using astronomy::KSPStockSystemFingerprint;
using astronomy::ParseTT;
using astronomy::StabilizeKSP;
using base::check_not_null;
using base::dynamic_cast_not_null;
using base::Error;
using base::FindOrDie;
using base::Fingerprint2011;
using base::HexadecimalEncoder;
using base::make_not_null_unique;
using base::not_null;
using base::OFStream;
using base::SerializeAsBytes;
using base::Status;
using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::DefinesFrame;
using geometry::EulerAngles;
using geometry::Identity;
using geometry::Normalize;
using geometry::Permutation;
using geometry::RigidTransformation;
using geometry::R3x3Matrix;
using geometry::Sign;
using physics::BarycentricRotatingDynamicFrame;
using physics::BodyCentredBodyDirectionDynamicFrame;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::BodySurfaceDynamicFrame;
using physics::BodySurfaceFrameField;
using physics::ComputeApsides;
using physics::ComputeNodes;
using physics::CoordinateFrameField;
using physics::DynamicFrame;
using physics::Frenet;
using physics::KeplerianElements;
using physics::MassiveBody;
using physics::RigidMotion;
using physics::SolarSystem;
using quantities::Force;
using quantities::Infinity;
using quantities::Length;
using quantities::MomentOfInertia;
using quantities::SIUnit;
using quantities::si::Kilogram;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using ::operator<<;

Plugin::Plugin(std::string const& game_epoch,
               std::string const& solar_system_epoch,
               Angle const& planetarium_rotation)
    : history_parameters_(DefaultHistoryParameters()),
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
  // parts, which have callbacks to remove themselves from |part_id_to_vessel_|,
  // which must therefore still exist.  This also causes the parts to be
  // destroyed, and therefore to destroy the pile-ups, which want to remove
  // themselves from |pile_up_|, which also exists.
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

void Plugin::InitializeEphemerisParameters(
    Ephemeris<Barycentric>::AccuracyParameters const& accuracy_parameters,
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters) {
  CHECK(initializing_);
  ephemeris_accuracy_parameters_ = accuracy_parameters;
  ephemeris_fixed_step_parameters_ = fixed_step_parameters;
}

void Plugin::InitializeHistoryParameters(
    Ephemeris<Barycentric>::FixedStepParameters const& parameters) {
  CHECK(initializing_);
  history_parameters_ = parameters;
}

void Plugin::InitializePsychohistoryParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const& parameters) {
  CHECK(initializing_);
  psychohistory_parameters_ = parameters;
}

void Plugin::EndInitialization() {
  CHECK(initializing_);
  SolarSystem<Barycentric> solar_system(gravity_model_, initial_state_);

  // If the system was constructed using keplerian elements, it may be the
  // stock KSP system in which case it needs to be stabilized.
  if (initial_state_.has_keplerian()) {
    auto const hierarchical_system = solar_system.MakeHierarchicalSystem();
    serialization::HierarchicalSystem message;
    hierarchical_system->WriteToMessage(&message);
    uint64_t const system_fingerprint =
        Fingerprint2011(SerializeAsBytes(message).get());
    LOG(INFO) << "System fingerprint is " << std::hex << std::uppercase
              << system_fingerprint;

    if (system_fingerprint == KSPStockSystemFingerprint) {
      LOG(WARNING) << "This appears to be the dreaded KSP stock system!";
      StabilizeKSP(solar_system);
      auto const hierarchical_system = solar_system.MakeHierarchicalSystem();
      serialization::HierarchicalSystem message;
      hierarchical_system->WriteToMessage(&message);
      uint64_t const system_fingerprint =
          Fingerprint2011(SerializeAsBytes(message).get());
      LOG(INFO) << "System fingerprint after stabilization is " << std::hex
                << std::uppercase << system_fingerprint;
      CHECK_EQ(KSPStabilizedSystemFingerprint, system_fingerprint)
          << "Attempt at stabilizing the KSP system failed!\n"
          << gravity_model_.DebugString() << "\n"
          << initial_state_.DebugString();
      LOG(INFO) << "This is the stabilized KSP system, all hail retrobop!";
    } else if (system_fingerprint == KSPStabilizedSystemFingerprint) {
      LOG(INFO) << "This is the stabilized KSP system, and we didn't have to "
                << "stabilize it ourselves.  All hail retrobop anyway!";
    } else {
      LOG(WARNING) << "This is an unknown system, we don't know anything about "
                   << "its stability:\n"
                   << gravity_model_.DebugString() << "\n"
                   << initial_state_.DebugString();
    }
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
  for (auto const& pair : celestials_) {
    Index const celestial_index = pair.first;
    auto const& celestial = pair.second;
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
  // TODO(egg): maybe these functions should take |Celestial*|s, and we should
  // then export |FindOrDie(celestials_, _)|.
  renderer_ = std::make_unique<Renderer>(
      sun_,
      make_not_null_unique<
          BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>>(
              ephemeris_.get(),
              sun_->body()));

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
  if (status.error() == Error::INVALID_ARGUMENT) {
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
  // |BodyWorld| with its y and z axes swapped (so that z is the polar axis).
  // The basis is right-handed.
  struct BodyFixed;
  Permutation<BodyWorld, BodyFixed> const body_mirror(
      Permutation<BodyWorld, BodyFixed>::XZY);

  auto const& body = *FindOrDie(celestials_, index)->body();

  OrthogonalMap<BodyWorld, World> const result =
      OrthogonalMap<WorldSun, World>::Identity() *
      sun_looking_glass.Inverse().Forget() *
      (PlanetariumRotation() *
       body.FromSurfaceFrame<BodyFixed>(current_time_)).Forget() *
      body_mirror.Forget();
  CHECK(result.Determinant().is_positive());
  return result.rotation();
}

Rotation<CelestialSphere, World> Plugin::CelestialSphereRotation()
    const {
  Permutation<CelestialSphere, Barycentric> const celestial_mirror(
      Permutation<CelestialSphere, Barycentric>::XZY);
  auto const result = OrthogonalMap<WorldSun, World>::Identity() *
                      sun_looking_glass.Inverse().Forget() *
                      PlanetariumRotation().Forget() *
                      celestial_mirror.Forget();
  CHECK(result.Determinant().is_positive());
  return result.rotation();
}

Angle Plugin::CelestialInitialRotation(Index const celestial_index) const {
  auto const& body = *FindOrDie(celestials_, celestial_index)->body();
  return body.AngleAt(game_epoch_);
}

Time Plugin::CelestialRotationPeriod(Index const celestial_index) const {
  auto const& body = *FindOrDie(celestials_, celestial_index)->body();
  // The result will be negative if the pole is the negative pole
  // (e.g. for Venus).  This is the convention KSP uses for retrograde rotation.
  return 2 * π * Radian / body.angular_frequency();
}

void Plugin::ClearWorldRotationalReferenceFrame() {
  angular_velocity_of_world_ = AngularVelocity<Barycentric>();
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
  auto it = vessels_.find(vessel_guid);
  if (it == vessels_.end()) {
    std::tie(it, inserted) =
        vessels_.emplace(
            vessel_guid,
            make_not_null_unique<Vessel>(vessel_guid,
                                         vessel_name,
                                         parent,
                                         ephemeris_.get(),
                                         DefaultPredictionParameters()));
  } else {
    inserted = false;
  }
  not_null<Vessel*> const vessel = it->second.get();
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
  RelativeDegreesOfFreedom<Barycentric> const relative =
      PlanetariumRotation().Inverse()(from_parent);
  ephemeris_->Prolong(current_time_);
  AddPart(
      vessel,
      part_id,
      name,
      InertiaTensor<RigidPart>::MakeWaterSphereInertiaTensor(1 * Kilogram),
      vessel->parent()->current_degrees_of_freedom(current_time_) + relative);
  // NOTE(egg): we do not keep the part; it may disappear just as we load, if
  // it happens to be a part with no physical significance (rb == null).
}

void Plugin::InsertOrKeepLoadedPart(
    PartId const part_id,
    std::string const& name,
    InertiaTensor<RigidPart> const& inertia_tensor,
    GUID const& vessel_guid,
    Index const main_body_index,
    DegreesOfFreedom<World> const& main_body_degrees_of_freedom,
    RigidMotion<RigidPart, World> const& part_rigid_motion,
    Time const& Δt) {
  not_null<Vessel*> const vessel = FindOrDie(vessels_, vessel_guid).get();
  CHECK(is_loaded(vessel));

  Instant const previous_time = current_time_ - Δt;
  OrthogonalMap<Barycentric, Barycentric> const Δplanetarium_rotation =
      Exp(Δt * angular_velocity_of_world_).Forget();
  // TODO(egg): Can we use |BarycentricToWorld| here?
  BodyCentredNonRotatingDynamicFrame<Barycentric, MainBodyCentred> const
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
  auto const world_to_barycentric_motion =
      main_body_frame.FromThisFrameAtTime(previous_time) *
      world_to_main_body_centred;

  auto it = part_id_to_vessel_.find(part_id);
  bool const part_found = it != part_id_to_vessel_.end();
  if (part_found) {
    not_null<Vessel*>& associated_vessel = it->second;
    not_null<Vessel*> const current_vessel = associated_vessel;
    if (vessel == current_vessel) {
    } else {
      associated_vessel = vessel;
      vessel->AddPart(current_vessel->ExtractPart(part_id));
    }
  } else {
    DegreesOfFreedom<World> const part_degrees_of_freedom(
        part_rigid_motion(DegreesOfFreedom<RigidPart>(RigidPart::origin,
                                                      Velocity<RigidPart>())));
    AddPart(vessel,
            part_id,
            name,
            inertia_tensor,
            world_to_barycentric_motion(part_degrees_of_freedom));
  }
  vessel->KeepPart(part_id);
  not_null<Part*> part = vessel->part(part_id);
  part->set_inertia_tensor(inertia_tensor);
}

void Plugin::IncrementPartIntrinsicForce(PartId const part_id,
                                         Vector<Force, World> const& force) {
  CHECK(!initializing_);
  not_null<Vessel*> const vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  vessel->part(part_id)->increment_intrinsic_force(
      renderer_->WorldToBarycentric(PlanetariumRotation())(force));
}

void Plugin::PrepareToReportCollisions() {
  for (auto const& [_, vessel] : vessels_) {
    // NOTE(egg): The lifetime requirement on the second argument of
    // |MakeSingleton| (which forwards to the argument of the constructor of
    // |Subset<Part>::Properties|) is that |part| outlives the constructed
    // |Properties|; since these are owned by |part|, this is true.
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
    if (kept_vessels_.erase(vessel)) {
      vessel->PrepareHistory(vessel_time);
      ++it;
    } else {
      loaded_vessels_.erase(vessel);
      LOG(INFO) << "Removing vessel " << vessel->ShortDebugString();
      renderer_->ClearTargetVesselIf(vessel);
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
    // vessel destroys its parts, which invalidates the intrusive |Subset| data
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
          history_parameters_,
          ephemeris_.get());
    });
  }
}

void Plugin::SetPartApparentDegreesOfFreedom(
    PartId const part_id,
    DegreesOfFreedom<World> const& degrees_of_freedom,
    DegreesOfFreedom<World> const& main_body_degrees_of_freedom) {
  // Define |ApparentBubble| as the reference frame with the axes of
  // |Barycentric| centred on the current main body.
  RigidMotion<World, ApparentBubble> world_to_apparent_bubble{
      RigidTransformation<World, ApparentBubble>{
          main_body_degrees_of_freedom.position(),
          ApparentBubble::origin,
          OrthogonalMap<Barycentric, ApparentBubble>::Identity() *
              renderer_->WorldToBarycentric(PlanetariumRotation())},
      renderer_->BarycentricToWorld(PlanetariumRotation())(
          -angular_velocity_of_world_),
      main_body_degrees_of_freedom.velocity()};

  not_null<Vessel*> vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  not_null<Part*> const part = vessel->part(part_id);
  CHECK(part->is_piled_up());
  part->containing_pile_up()->SetPartApparentDegreesOfFreedom(
      part, world_to_apparent_bubble(degrees_of_freedom));
}

DegreesOfFreedom<World> Plugin::GetPartActualDegreesOfFreedom(
    PartId const part_id,
    RigidMotion<Barycentric, World> const& barycentric_to_world) const {
  return barycentric_to_world(
             FindOrDie(part_id_to_vessel_,
                       part_id)->part(part_id)->degrees_of_freedom());
}

DegreesOfFreedom<World> Plugin::CelestialWorldDegreesOfFreedom(
    Index const index,
    RigidMotion<Barycentric, World> const& barycentric_to_world,
    Instant const& time) const {
  return barycentric_to_world(
             FindOrDie(celestials_, index)->
                 trajectory().EvaluateDegreesOfFreedom(time));
}

RigidMotion<Barycentric, World> Plugin::BarycentricToWorld(
    bool const reference_part_is_unmoving,
    PartId const reference_part_id,
    std::optional<Position<World>> const& main_body_centre) const {
  BodyCentredNonRotatingDynamicFrame<Barycentric, MainBodyCentred> const
      main_body_frame{ephemeris_.get(), main_body_};

  auto const barycentric_to_main_body_motion =
      main_body_frame.ToThisFrameAtTime(current_time_);
  // In coordinates, this rotation is the identity.
  auto const barycentric_to_main_body_rotation =
      barycentric_to_main_body_motion.rigid_transformation().linear_map();
  auto const reference_part_degrees_of_freedom =
      barycentric_to_main_body_motion(
          FindOrDie(part_id_to_vessel_, reference_part_id)->
              part(reference_part_id)->degrees_of_freedom());

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
      /*velocity_of_to_frame_origin=*/Velocity<World>()}.Inverse();
    }
  }();
  return main_body_to_world * barycentric_to_main_body_motion;
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  CHECK(!initializing_);
  CHECK_GT(t, current_time_);

  for (not_null<Vessel*> const vessel : loaded_vessels_) {
    vessel->ClearAllIntrinsicForces();
  }

  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
  ephemeris_->Prolong(current_time_);
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
    if (vessel->psychohistory().back().time < current_time_) {
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
        Status const status = pile_up->DeformAndAdvanceTime(current_time_);
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
  Status const status = future.get();
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

void Plugin::ForgetAllHistoriesBefore(Instant const& t) const {
  CHECK(!initializing_);
  CHECK_LT(t, current_time_);
  ephemeris_->EventuallyForgetBefore(t);
  for (auto const& [_, vessel] : vessels_) {
    vessel->ForgetBefore(t);
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
      vessel->psychohistory().back().degrees_of_freedom -
      vessel->parent()->current_degrees_of_freedom(current_time_);
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  return result;
}

RelativeDegreesOfFreedom<AliceSun> Plugin::CelestialFromParent(
    Index const celestial_index) const {
  CHECK(!initializing_);
  ephemeris_->Prolong(current_time_);
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
  if (renderer_->HasTargetVessel()) {
    renderer_->GetTargetVessel().set_prediction_adaptive_step_parameters(
        prediction_adaptive_step_parameters);
  }
  FindOrDie(vessels_, vessel_guid)
      ->set_prediction_adaptive_step_parameters(
          prediction_adaptive_step_parameters);
}

void Plugin::UpdatePrediction(GUID const& vessel_guid) const {
  CHECK(!initializing_);
  Vessel& vessel = *FindOrDie(vessels_, vessel_guid);

  // If there is a target vessel, ensure that the prediction of |vessel| is not
  // longer than that of the target vessel.  This is necessary to build the
  // targetting frame.
  if (renderer_->HasTargetVessel()) {
    Vessel& target_vessel = renderer_->GetTargetVessel();
    target_vessel.RefreshPrediction();
    vessel.RefreshPrediction(target_vessel.prediction().back().time);
  } else {
    vessel.RefreshPrediction();
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

void Plugin::ComputeAndRenderApsides(
    Index const celestial_index,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    int const max_points,
    std::unique_ptr<DiscreteTrajectory<World>>& apoapsides,
    std::unique_ptr<DiscreteTrajectory<World>>& periapsides) const {
  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(FindOrDie(celestials_, celestial_index)->trajectory(),
                 begin,
                 end,
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

void Plugin::ComputeAndRenderClosestApproaches(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    int const max_points,
    std::unique_ptr<DiscreteTrajectory<World>>& closest_approaches) const {
  CHECK(renderer_->HasTargetVessel());

  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(renderer_->GetTargetVessel().prediction(),
                 begin,
                 end,
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
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    int const max_points,
    std::unique_ptr<DiscreteTrajectory<World>>& ascending,
    std::unique_ptr<DiscreteTrajectory<World>>& descending) const {
  auto const trajectory_in_plotting =
      renderer_->RenderBarycentricTrajectoryInPlotting(begin, end);

  auto const* const cast_plotting_frame = dynamic_cast<
      BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation> const*>(
      &*renderer_->GetPlottingFrame());
  // The body-centred non rotating frame does not rotate; its reference plane is
  // not an inherent dynamical property.  When using this frame, discard the
  // nodes if they are far enough from the central body that its equator is
  // irrelevant.
  Length const threshold =
      cast_plotting_frame == nullptr
          ? Infinity<Length>()
          : EquatorRelevanceThreshold(
                *dynamic_cast_not_null<RotatingBody<Barycentric> const*>(
                    cast_plotting_frame->centre()));
  auto const show_node = [threshold](DegreesOfFreedom<Navigation> const& dof) {
    return (dof.position() - Navigation::origin).Norm() < threshold;
  };

  DiscreteTrajectory<Navigation> ascending_trajectory;
  DiscreteTrajectory<Navigation> descending_trajectory;
  // The so-called North is orthogonal to the plane of the trajectory.
  ComputeNodes(trajectory_in_plotting->begin(),
               trajectory_in_plotting->end(),
               Vector<double, Navigation>({0, 0, 1}),
               max_points,
               ascending_trajectory,
               descending_trajectory,
               show_node);

  ascending = renderer_->RenderPlottingTrajectoryInWorld(
                  current_time_,
                  ascending_trajectory.begin(),
                  ascending_trajectory.end(),
                  sun_world_position,
                  PlanetariumRotation());
  descending = renderer_->RenderPlottingTrajectoryInWorld(
                   current_time_,
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

not_null<std::unique_ptr<Planetarium>> Plugin::NewPlanetarium(
    Planetarium::Parameters const& parameters,
    Perspective<Navigation, Camera> const& perspective)
    const {
  return make_not_null_unique<Planetarium>(parameters,
                                           perspective,
                                           ephemeris_.get(),
                                           renderer_->GetPlottingFrame());
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBarycentricRotatingNavigationFrame(
    Index const primary_index,
    Index const secondary_index) const {
  CHECK(!initializing_);
  // TODO(egg): these should be const, use a custom comparator in the map.
  Celestial const& primary = *FindOrDie(celestials_, primary_index);
  Celestial const& secondary = *FindOrDie(celestials_, secondary_index);
  return make_not_null_unique<
      BarycentricRotatingDynamicFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          primary.body(),
          secondary.body());
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBodyCentredBodyDirectionNavigationFrame(
    Index const primary_index,
    Index const secondary_index) const {
  CHECK(!initializing_);
  // TODO(egg): these should be const, use a custom comparator in the map.
  Celestial const& primary = *FindOrDie(celestials_, primary_index);
  Celestial const& secondary = *FindOrDie(celestials_, secondary_index);
  return make_not_null_unique<
      BodyCentredBodyDirectionDynamicFrame<Barycentric, Navigation>>(
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
      BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          reference_body.body());
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBodySurfaceNavigationFrame(
    Index const reference_body_index) const {
  CHECK(!initializing_);
  Celestial const& reference_body =
      *FindOrDie(celestials_, reference_body_index);
  return make_not_null_unique<BodySurfaceDynamicFrame<Barycentric, Navigation>>(
      ephemeris_.get(),
      reference_body.body());
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

  struct RightHandedNavball;

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
      plugin_->ephemeris_->Prolong(current_time);

      Position<Navigation> const q_in_plotting =
          renderer.WorldToPlotting(current_time,
                                   sun_world_position_,
                                   planetarium_rotation)(q);

      OrthogonalMap<RightHandedNavball, Barycentric> const
          right_handed_navball_to_barycentric =
              barycentric_right_handed_field_ == nullptr
                  ? renderer.PlottingToBarycentric(current_time) *
                        navigation_right_handed_field_->
                            FromThisFrame(q_in_plotting).Forget()
                  : barycentric_right_handed_field_->FromThisFrame(
                        renderer.WorldToBarycentric(
                            current_time,
                            sun_world_position_,
                            planetarium_rotation)(q)).Forget();

      // KSP's navball has x west, y up, z south.
      // We want x north, y east, z down.
      OrthogonalMap<Navball, World> const orthogonal_map =
          renderer.BarycentricToWorld(planetarium_rotation) *
          right_handed_navball_to_barycentric *
          Permutation<World, RightHandedNavball>(
              Permutation<World, RightHandedNavball>::XZY).Forget() *
          Rotation<Navball, World>(π / 2 * Radian,
                                   Bivector<double, World>({0, 1, 0}),
                                   DefinesFrame<Navball>()).Forget();
      CHECK(orthogonal_map.Determinant().is_positive());
      return orthogonal_map.rotation();
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
  auto* const plotting_frame_as_body_surface_dynamic_frame =
      dynamic_cast<BodySurfaceDynamicFrame<Barycentric, Navigation> const*>(
          &*renderer_->GetPlottingFrame());
  if (plotting_frame_as_body_surface_dynamic_frame == nullptr) {
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
                plotting_frame_as_body_surface_dynamic_frame->centre()),
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
  auto const back = vessel.psychohistory().back();
  return VesselVelocity(back.time, back.degrees_of_freedom);
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

Renderer& Plugin::renderer() {
  return *renderer_;
}

Renderer const& Plugin::renderer() const {
  return *renderer_;
}

void Plugin::WriteToMessage(
    not_null<serialization::Plugin*> const message) const {
  LOG(INFO) << __FUNCTION__;
  CHECK(!initializing_);
  ephemeris_->Prolong(current_time_);
  std::map<not_null<Celestial const*>, Index const> celestial_to_index;
  for (auto const& pair : celestials_) {
    Index const index = pair.first;
    auto const& owned_celestial = pair.second;
    celestial_to_index.emplace(owned_celestial.get(), index);
  }
  for (auto const& pair : celestials_) {
    Index const index = pair.first;
    auto const& owned_celestial = pair.second.get();
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
  for (auto const& pair : part_id_to_vessel_) {
    PartId const part_id = pair.first;
    not_null<Vessel*> const vessel = pair.second;
    (*message->mutable_part_id_to_vessel())[part_id] = vessel_to_guid[vessel];
  }

  ephemeris_->WriteToMessage(message->mutable_ephemeris());

  history_parameters_.WriteToMessage(message->mutable_history_parameters());
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

  plugin->game_epoch_ = Instant::ReadFromMessage(message.game_epoch());
  plugin->current_time_ = Instant::ReadFromMessage(message.current_time());
  plugin->planetarium_rotation_ =
      Angle::ReadFromMessage(message.planetarium_rotation());

  // The ephemeris constructed here is *not* prolonged and needs to be
  // explicitly prolonged to cover all the instants that we care about.
  plugin->ephemeris_ =
      Ephemeris<Barycentric>::ReadFromMessage(message.ephemeris());
  plugin->ephemeris_->Prolong(plugin->game_epoch_);
  plugin->ephemeris_->Prolong(plugin->current_time_);

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

  for (auto const& pair : message.part_id_to_vessel()) {
    PartId const part_id = pair.first;
    GUID const guid = pair.second;
    auto const& vessel = FindOrDie(plugin->vessels_, guid);
    plugin->part_id_to_vessel_.emplace(part_id, vessel.get());
  }

  plugin->sun_ = FindOrDie(plugin->celestials_, message.sun_index()).get();
  plugin->main_body_ = plugin->sun_->body();
  plugin->UpdatePlanetariumRotation();

  bool const is_pre_cauchy = message.has_pre_cauchy_plotting_frame();
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
    Ephemeris<Barycentric>::FixedStepParameters const& history_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        psychohistory_parameters)
    : history_parameters_(history_parameters),
      psychohistory_parameters_(psychohistory_parameters),
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
  // The z axis of |PlanetariumFrame| is the pole of |main_body_|, and its x
  // axis is the origin of body rotation (the intersection between the
  // |Barycentric| xy plane and the plane of |main_body_|'s equator, or the y
  // axis of |Barycentric| if they coincide).
  // This can be expressed using Euler angles, see figures 1 and 2 of
  // http://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2009reprint.pdf.
  struct PlanetariumFrame;

  CHECK_NOTNULL(main_body_);
  Rotation<Barycentric, PlanetariumFrame> const to_planetarium(
      π / 2 * Radian + main_body_->right_ascension_of_pole(),
      π / 2 * Radian - main_body_->declination_of_pole(),
      0 * Radian,
      EulerAngles::ZXZ,
      DefinesFrame<PlanetariumFrame>{});
  cached_planetarium_rotation_ =
      Rotation<PlanetariumFrame, AliceSun>(
          planetarium_rotation_,
          Bivector<double, PlanetariumFrame>({0, 0, 1}),
          DefinesFrame<AliceSun>{}) *
      to_planetarium;
}

Velocity<World> Plugin::VesselVelocity(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) const {
  DegreesOfFreedom<Navigation> const plotting_frame_degrees_of_freedom =
      renderer_->BarycentricToPlotting(time)(degrees_of_freedom);
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

void Plugin::AddPart(not_null<Vessel*> const vessel,
                     PartId const part_id,
                     std::string const& name,
                     InertiaTensor<RigidPart> const& inertia_tensor,
                     DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  auto const [it, inserted] = part_id_to_vessel_.emplace(part_id, vessel);
  CHECK(inserted) << NAMED(part_id);
  auto deletion_callback = [it = it, &map = part_id_to_vessel_] {
    map.erase(it);
  };
  auto part = make_not_null_unique<Part>(part_id,
                                         name,
                                         inertia_tensor,
                                         degrees_of_freedom,
                                         std::move(deletion_callback));
  vessel->AddPart(std::move(part));
}

bool Plugin::is_loaded(not_null<Vessel*> vessel) const {
  return Contains(loaded_vessels_, vessel);
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
