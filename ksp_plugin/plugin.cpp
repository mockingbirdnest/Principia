
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <experimental/filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <set>

#include "base/file.hpp"
#include "base/hexadecimal.hpp"
#include "base/map_util.hpp"
#include "base/not_null.hpp"
#include "base/optional_logging.hpp"
#include "base/unique_ptr_logging.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/identity.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/permutation.hpp"
#include "glog/logging.h"
#include "glog/stl_logging.h"
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

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using base::dynamic_cast_not_null;
using base::Error;
using base::FindOrDie;
using base::Fingerprint2011;
using base::FingerprintCat2011;
using base::make_not_null_unique;
using base::OFStream;
using base::not_null;
using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::DefinesFrame;
using geometry::EulerAngles;
using geometry::Identity;
using geometry::Normalize;
using geometry::Permutation;
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
using physics::RigidTransformation;
using quantities::Force;
using quantities::Length;
using quantities::si::Kilogram;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using ::operator<<;

namespace {

Length const fitting_tolerance = 1 * Milli(Metre);

std::uint64_t const ksp_stock_system_fingerprint = 0x025779971BA2BFD7u;
std::uint64_t const ksp_fixed_system_fingerprint = 0x1248ADFCBD8BCE64u;

// The map between the vector spaces of |WorldSun| and |AliceSun|.
Permutation<WorldSun, AliceSun> const sun_looking_glass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

}  // namespace

Plugin::Plugin(Instant const& game_epoch,
               Instant const& solar_system_epoch,
               Angle const& planetarium_rotation)
    : history_parameters_(DefaultHistoryParameters()),
      prolongation_parameters_(DefaultProlongationParameters()),
      prediction_parameters_(DefaultPredictionParameters()),
      planetarium_rotation_(planetarium_rotation),
      game_epoch_(game_epoch),
      current_time_(solar_system_epoch) {}

void Plugin::InsertCelestialAbsoluteCartesian(
    Index const celestial_index,
    std::experimental::optional<Index> const& parent_index,
    DegreesOfFreedom<Barycentric> const& initial_state,
    not_null<std::unique_ptr<RotatingBody<Barycentric> const>> body) {
  LOG(INFO) << __FUNCTION__ << "\n"
            << NAMED(celestial_index) << "\n"
            << NAMED(parent_index) << "\n"
            << NAMED(initial_state) << "\n"
            << NAMED(body);
  CHECK(initializing_) << "Celestial bodies should be inserted before the end "
                       << "of initialization";
  CHECK(!hierarchical_initialization_);
  if (!absolute_initialization_) {
    absolute_initialization_.emplace();
  }
  auto const inserted = celestials_.emplace(
      celestial_index,
      std::make_unique<Celestial>(body.get()));
  CHECK(inserted.second) << "Body already exists at index " << celestial_index;
  not_null<Celestial*> const celestial = inserted.first->second.get();
  absolute_initialization_->bodies.emplace(celestial_index, std::move(body));
  if (parent_index) {
    not_null<Celestial const*> parent =
        FindOrDie(celestials_, *parent_index).get();
    celestial->set_parent(parent);
  } else {
    CHECK(sun_ == nullptr);
    sun_ = celestial;
  }
  absolute_initialization_->initial_state.emplace(celestial_index,
                                                  initial_state);
}

Plugin::~Plugin() {
  // We must manually destroy the vessels, triggering the destruction of the
  // parts, which have callbacks to remove themselves from |part_id_to_vessel_|,
  // which must therefore still exist.  This also removes the parts from the
  // pile-ups, which also exist.
  vessels_.clear();
}

void Plugin::InsertCelestialJacobiKeplerian(
    Index const celestial_index,
    std::experimental::optional<Index> const& parent_index,
    std::experimental::optional<KeplerianElements<Barycentric>> const&
        keplerian_elements,
    not_null<std::unique_ptr<RotatingBody<Barycentric>>> body) {
  LOG(INFO) << __FUNCTION__ << "\n"
            << NAMED(celestial_index) << "\n"
            << NAMED(parent_index) << "\n"
            << NAMED(keplerian_elements) << "\n"
            << NAMED(body);
  CHECK(initializing_) << "Celestial bodies should be inserted before the end "
                       << "of initialization";
  CHECK(!absolute_initialization_);
  CHECK_EQ((bool)parent_index, (bool)keplerian_elements);
  CHECK_EQ((bool)parent_index, (bool)hierarchical_initialization_);
  RotatingBody<Barycentric>* const unowned_body = body.get();
  if (hierarchical_initialization_) {
    hierarchical_initialization_->system.Add(
        std::move(body),
        hierarchical_initialization_->indices_to_bodies[*parent_index],
        *keplerian_elements);
  } else {
    hierarchical_initialization_.emplace(std::move(body));
  }
  bool inserted =
      hierarchical_initialization_->parents.emplace(celestial_index,
                                                    parent_index).second;
  inserted &=
      hierarchical_initialization_->
          indices_to_bodies.emplace(celestial_index, unowned_body).second;
  CHECK(inserted);

  // Record the fingerprints of the parameters to detect if we are in KSP stock.
  CHECK(celestial_jacobi_keplerian_fingerprints_.insert(
            FingerprintCelestialJacobiKeplerian(celestial_index,
                                                parent_index,
                                                keplerian_elements,
                                                *unowned_body)).second);
}

void Plugin::EndInitialization() {
  CHECK(initializing_);
  if (hierarchical_initialization_) {
    std::uint64_t system_fingerprint = 0;
    for (std::uint64_t fingerprint : celestial_jacobi_keplerian_fingerprints_) {
      system_fingerprint = FingerprintCat2011(system_fingerprint, fingerprint);
    }
    LOG(INFO) << "System fingerprint is " << std::hex << system_fingerprint;
    if (system_fingerprint == ksp_stock_system_fingerprint) {
      is_ksp_stock_system_ = true;
      LOG(WARNING) << "This appears to be the dreaded KSP stock system!";
    } else if (system_fingerprint == ksp_fixed_system_fingerprint) {
      LOG(INFO) << "This is the fixed KSP system, all hail retrobop!";
    }

    HierarchicalSystem<Barycentric>::BarycentricSystem system =
        hierarchical_initialization_->system.ConsumeBarycentricSystem();
    std::map<not_null<RotatingBody<Barycentric> const*>,
             Index> bodies_to_indices;
    for (auto const& index_body :
             hierarchical_initialization_->indices_to_bodies) {
      bodies_to_indices[index_body.second] = index_body.first;
    }
    auto const parents = std::move(hierarchical_initialization_->parents);
    hierarchical_initialization_ = std::experimental::nullopt;
#if LOG_KSP_SYSTEM
    OFStream file;
    if (system_fingerprint == ksp_stock_system_fingerprint) {
      file = OFStream(TEMP_DIR / "ksp_stock_system.proto.hex");
    } else if (system_fingerprint == ksp_fixed_system_fingerprint) {
      file = OFStream(TEMP_DIR / "ksp_fixed_system.proto.hex");
    } else {
      file = OFStream(TEMP_DIR / "unknown_system.proto.hex");
    }
    std::string bytes;
    base::UniqueArray<std::uint8_t> hex;
#endif
    for (int i = 0; i < system.bodies.size(); ++i) {
#if LOG_KSP_SYSTEM
      serialization::MassiveBody body_message;
      serialization::Pair degrees_of_freedom_message;
      system.bodies[i]->WriteToMessage(&body_message);
      system.degrees_of_freedom[i].WriteToMessage(&degrees_of_freedom_message);
      body_message.SerializeToString(&bytes);
      hex = base::UniqueArray<std::uint8_t>((bytes.size() << 1) + 1);
      base::HexadecimalEncode(
          base::Array<std::uint8_t const>(
              reinterpret_cast<std::uint8_t const*>(bytes.data()),
              bytes.size()),
          hex.get());
      hex.data[hex.size - 1] = 0;
      file << reinterpret_cast<char const*>(hex.data.get()) << "\n";
      degrees_of_freedom_message.SerializeToString(&bytes);
      hex = base::UniqueArray<std::uint8_t>((bytes.size() << 1) + 1);
      base::HexadecimalEncode(
          base::Array<std::uint8_t const>(
              reinterpret_cast<std::uint8_t const*>(bytes.data()),
              bytes.size()),
          hex.get());
      hex.data[hex.size - 1] = 0;
      file << reinterpret_cast<char const*>(hex.data.get()) << "\n";
#endif
      auto rotating_body = dynamic_cast_not_null<
          std::unique_ptr<RotatingBody<Barycentric> const>>(
              std::move(system.bodies[i]));
      Index const celestial_index = bodies_to_indices[rotating_body.get()];
      InsertCelestialAbsoluteCartesian(
          celestial_index,
          FindOrDie(parents, celestial_index),
          system.degrees_of_freedom[i],
          std::move(rotating_body));
    }
  }
  CHECK(absolute_initialization_);
  CHECK_NOTNULL(sun_);
  main_body_ = sun_->body();
  UpdatePlanetariumRotation();
  initializing_.Flop();

  InitializeEphemerisAndSetCelestialTrajectories();

  // Log the serialized ephemeris.
  serialization::Ephemeris ephemeris_message;
  ephemeris_->WriteToMessage(&ephemeris_message);
  std::string const bytes = ephemeris_message.SerializeAsString();
  base::UniqueArray<std::uint8_t> const hex((bytes.size() << 1) + 1);
  base::HexadecimalEncode(
      base::Array<std::uint8_t const>(
          reinterpret_cast<std::uint8_t const*>(bytes.data()), bytes.size()),
      hex.get());
  hex.data[hex.size - 1] = 0;
  // Begin and end markers to make sure the hex did not get clipped (this might
  // happen if the message is very big).
  LOG(INFO) << "Ephemeris at initialization:\nbegin\n"
            << reinterpret_cast<char const*>(hex.data.get()) << "\nend";
}

bool Plugin::IsKspStockSystem() const {
  CHECK(!initializing_);
  return is_ksp_stock_system_;
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
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(celestial_index) << '\n' << NAMED(parent_index);
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
  CHECK(result.Determinant().Positive());
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
  CHECK(result.Determinant().Positive());
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

void Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                std::string const& vessel_name,
                                Index const parent_index,
                                bool const loaded,
                                bool& inserted) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(vessel_guid) << '\n'
          << NAMED(parent_index) << '\n'
          << NAMED(loaded) << '\n';
  CHECK(!initializing_);
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  auto it = vessels_.find(vessel_guid);
  if (it == vessels_.end()) {
    std::tie(it, inserted) =
        vessels_.emplace(vessel_guid,
                         make_not_null_unique<Vessel>(vessel_guid,
                                                      vessel_name,
                                                      parent,
                                                      ephemeris_.get(),
                                                      prediction_parameters_));
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
  not_null<Vessel*> const vessel =
      find_vessel_by_guid_or_die(vessel_guid).get();
  RelativeDegreesOfFreedom<Barycentric> const relative =
      PlanetariumRotation().Inverse()(from_parent);
  ephemeris_->Prolong(current_time_);
  AddPart(vessel,
          part_id,
          name,
          1 * Kilogram,
          vessel->parent()->current_degrees_of_freedom(current_time_) +
              relative);
  // NOTE(egg): we do not keep the part; it may disappear just as we load, if
  // it happens to be a part with no physical significance (rb == null).
}

void Plugin::InsertOrKeepLoadedPart(
    PartId const part_id,
    std::string const& name,
    Mass const& mass,
    GUID const& vessel_guid,
    Index const main_body_index,
    DegreesOfFreedom<World> const& main_body_degrees_of_freedom,
    DegreesOfFreedom<World> const& part_degrees_of_freedom) {
  not_null<Vessel*> const vessel =
      find_vessel_by_guid_or_die(vessel_guid).get();
  CHECK(is_loaded(vessel));

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
    enum class LocalTag { tag };
    using MainBodyCentred =
        geometry::Frame<LocalTag, LocalTag::tag, /*frame_is_inertial=*/false>;
    BodyCentredNonRotatingDynamicFrame<Barycentric, MainBodyCentred> const
        main_body_frame{ephemeris_.get(),
                        FindOrDie(celestials_, main_body_index)->body()};
    RigidMotion<World, MainBodyCentred> const world_to_main_body_centred{
        RigidTransformation<World, MainBodyCentred>{
            main_body_degrees_of_freedom.position(),
            MainBodyCentred::origin,
            main_body_frame.ToThisFrameAtTime(current_time_).orthogonal_map() *
                WorldToBarycentric()},
        AngularVelocity<World>(),
        main_body_degrees_of_freedom.velocity()};
    auto const world_to_barycentric =
        main_body_frame.FromThisFrameAtTime(current_time_) *
        world_to_main_body_centred;

    AddPart(vessel,
            part_id,
            name,
            mass,
            world_to_barycentric(part_degrees_of_freedom));
  }
  vessel->KeepPart(part_id);
  not_null<Part*> part = vessel->part(part_id);
  part->set_mass(mass);
}

void Plugin::IncrementPartIntrinsicForce(PartId const part_id,
                                         Vector<Force, World> const& force) {
  CHECK(!initializing_);
  not_null<Vessel*> const vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  vessel->part(part_id)->increment_intrinsic_force(WorldToBarycentric()(force));
}

void Plugin::PrepareToReportCollisions() {
  for (auto const& pair : vessels_) {
    Vessel& vessel = *pair.second;
    // TODO(egg): we're taking the address of a parameter passed by reference
    // here; but then I don't think I want to pass this by pointer, it's quite
    // convenient everywhere else...
    vessel.ForAllParts(
        [](Part& part) { Subset<Part>::MakeSingleton(part, &part); });
  }
}

void Plugin::ReportCollision(PartId const part1, PartId const part2) const {
  Part& p1 = *FindOrDie(part_id_to_vessel_, part1)->part(part1);
  Part& p2 = *FindOrDie(part_id_to_vessel_, part2)->part(part2);
  Subset<Part>::Unite(Subset<Part>::Find(p1), Subset<Part>::Find(p2));
}

void Plugin::FreeVesselsAndPartsAndCollectPileUps() {
  CHECK(!initializing_);

  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    not_null<Vessel*> vessel = it->second.get();
    if (kept_vessels_.erase(vessel)) {
      vessel->PreparePsychohistory(current_time_);
      ++it;
    } else {
      CHECK(!is_loaded(vessel));
      LOG(INFO) << "Removing vessel " << vessel->ShortDebugString();
      if (target_ && target_->vessel == vessel) {
        target_ = std::experimental::nullopt;
      }
      it = vessels_.erase(it);
    }
  }
  CHECK(kept_vessels_.empty());

  // Free old parts.
  for (not_null<Vessel*> const vessel : loaded_vessels_) {
    vessel->FreeParts();
  }

  // Bind the vessels.
  for (auto const& pair : vessels_) {
    Vessel& vessel = *pair.second;
    vessel.ForSomePart([&vessel, this](Part& first_part) {
      vessel.ForAllParts([&first_part](Part& part) {
        Subset<Part>::Unite(Subset<Part>::Find(first_part),
                            Subset<Part>::Find(part));
      });
    });
  }

  // We only need to collect one part per vessel, since the other parts are in
  // the same subset.
  for (auto const& pair : vessels_) {
    Vessel& vessel = *pair.second;
    vessel.ForSomePart([this](Part& first_part) {
      Subset<Part>::Find(first_part).mutable_properties().Collect(
          &pile_ups_,
          current_time_,
          DefaultProlongationParameters(),
          DefaultHistoryParameters(),
          ephemeris_.get());
    });
  }
}

void Plugin::SetPartApparentDegreesOfFreedom(
    PartId const part_id,
    DegreesOfFreedom<World> const& degrees_of_freedom) {
  RigidMotion<World, ApparentBubble> world_to_apparent_bubble{
      RigidTransformation<World, ApparentBubble>{
          World::origin,
          ApparentBubble::origin,
          OrthogonalMap<Barycentric, ApparentBubble>::Identity() *
              WorldToBarycentric()},
      AngularVelocity<World>{},
      Velocity<World>{}};
  not_null<Vessel*> vessel = FindOrDie(part_id_to_vessel_, part_id);
  CHECK(is_loaded(vessel));
  not_null<Part*> const part = vessel->part(part_id);
  CHECK(part->is_piled_up());
  part->containing_pile_up()->iterator()->SetPartApparentDegreesOfFreedom(
      part, world_to_apparent_bubble(degrees_of_freedom));
}

void Plugin::AdvanceParts(Instant const& t) {
  CHECK(!initializing_);
  CHECK_GT(t, current_time_);

  ephemeris_->Prolong(t);
  for (PileUp& pile_up : pile_ups_) {
    pile_up.DeformPileUpIfNeeded();
    pile_up.AdvanceTime(t);
    // TODO(egg): now that |NudgeParts| doesn't need the bubble barycentre
    // anymore, it could be part of |PileUp::AdvanceTime|.
    pile_up.NudgeParts();
  }
}

DegreesOfFreedom<World> Plugin::GetPartActualDegreesOfFreedom(
    PartId const part_id,
    PartId const part_at_origin) const {
  auto const world_origin = FindOrDie(part_id_to_vessel_, part_at_origin)->
                                part(part_at_origin)->
                                degrees_of_freedom();
  RigidMotion<Barycentric, World> barycentric_to_world{
      RigidTransformation<Barycentric, World>{
          world_origin.position(), World::origin, BarycentricToWorld()},
      AngularVelocity<Barycentric>{},
      world_origin.velocity()};
  return barycentric_to_world(FindOrDie(part_id_to_vessel_, part_id)->
                                  part(part_id)->
                                  degrees_of_freedom());
}

DegreesOfFreedom<World> Plugin::CelestialWorldDegreesOfFreedom(
    Index const index,
    PartId const part_at_origin) const {
  auto const world_origin = FindOrDie(part_id_to_vessel_, part_at_origin)->
                                part(part_at_origin)->
                                degrees_of_freedom();
  RigidMotion<Barycentric, World> barycentric_to_world{
      RigidTransformation<Barycentric, World>{
          world_origin.position(), World::origin, BarycentricToWorld()},
      AngularVelocity<Barycentric>{},
      world_origin.velocity()};
  return barycentric_to_world(
      FindOrDie(celestials_, index)->
          trajectory().EvaluateDegreesOfFreedom(current_time_));
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(t) << '\n' << NAMED(planetarium_rotation);
  CHECK(!initializing_);
  CHECK_GT(t, current_time_);

  if (!vessels_.empty()) {
    bool tails_are_empty;
    vessels_.begin()->second->ForSomePart([&tails_are_empty](Part& part) {
      tails_are_empty = part.tail().Empty();
    });
    if (tails_are_empty) {
      AdvanceParts(t);
    }
  }

  for (auto const& pair : vessels_) {
    Vessel& vessel = *pair.second;
    vessel.AdvanceTime();
  }
  for (not_null<Vessel*> const vessel : loaded_vessels_) {
    vessel->ClearAllIntrinsicForces();
  }

  VLOG(1) << "Time has been advanced" << '\n'
          << "from : " << current_time_ << '\n'
          << "to   : " << t;
  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
  UpdatePlanetariumRotation();
  loaded_vessels_.clear();
}

void Plugin::ForgetAllHistoriesBefore(Instant const& t) const {
  CHECK(!initializing_);
  CHECK_LT(t, current_time_);
  ephemeris_->ForgetBefore(t);
  for (auto const& pair : vessels_) {
    not_null<std::unique_ptr<Vessel>> const& vessel = pair.second;
    vessel->ForgetBefore(t);
  }
}

RelativeDegreesOfFreedom<AliceSun> Plugin::VesselFromParent(
    Index const parent_index,
    GUID const& vessel_guid) const {
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  if (vessel->parent() != parent) {
    vessel->set_parent(parent);
  }
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      vessel->psychohistory().last().degrees_of_freedom() -
      vessel->parent()->current_degrees_of_freedom(current_time_);
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  VLOG(1) << "Vessel " << vessel->ShortDebugString()
          << " is at parent degrees of freedom + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
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
  VLOG(1) << "Celestial at index " << celestial_index
          << " is at parent degrees of freedom + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

void Plugin::UpdatePrediction(GUID const& vessel_guid) const {
  CHECK(!initializing_);
  find_vessel_by_guid_or_die(vessel_guid)->UpdatePrediction(
      current_time_ + prediction_length_);
}

void Plugin::CreateFlightPlan(GUID const& vessel_guid,
                              Instant const& final_time,
                              Mass const& initial_mass) const {
  CHECK(!initializing_);
  find_vessel_by_guid_or_die(vessel_guid)->CreateFlightPlan(
      final_time,
      initial_mass,
      prediction_parameters_);
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
Plugin::RenderBarycentricTrajectoryInWorld(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position) const {
  auto const trajectory_in_navigation =
      RenderBarycentricTrajectoryInNavigation(begin, end);
  auto trajectory_in_world =
      RenderNavigationTrajectoryInWorld(trajectory_in_navigation->Begin(),
                                        trajectory_in_navigation->End(),
                                        sun_world_position);
  return trajectory_in_world;
}

void Plugin::ComputeAndRenderApsides(
    Index const celestial_index,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    std::unique_ptr<DiscreteTrajectory<World>>& apoapsides,
    std::unique_ptr<DiscreteTrajectory<World>>& periapsides) const {
  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(FindOrDie(celestials_, celestial_index)->trajectory(),
                 begin,
                 end,
                 apoapsides_trajectory,
                 periapsides_trajectory);
  apoapsides =
      RenderBarycentricTrajectoryInWorld(apoapsides_trajectory.Begin(),
                                         apoapsides_trajectory.End(),
                                         sun_world_position);
  periapsides =
      RenderBarycentricTrajectoryInWorld(periapsides_trajectory.Begin(),
                                         periapsides_trajectory.End(),
                                         sun_world_position);
}

void Plugin::ComputeAndRenderClosestApproaches(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    std::unique_ptr<DiscreteTrajectory<World>>& closest_approaches) const {
  CHECK(target_);

  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(target_->vessel->prediction(),
                 begin,
                 end,
                 apoapsides_trajectory,
                 periapsides_trajectory);
  closest_approaches =
      RenderBarycentricTrajectoryInWorld(periapsides_trajectory.Begin(),
                                         periapsides_trajectory.End(),
                                         sun_world_position);
}

void Plugin::ComputeAndRenderNodes(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    std::unique_ptr<DiscreteTrajectory<World>>& ascending,
    std::unique_ptr<DiscreteTrajectory<World>>& descending) const {
  CHECK(target_);
  auto const trajectory_in_navigation =
      RenderBarycentricTrajectoryInNavigation(begin, end);
  DiscreteTrajectory<Navigation> ascending_trajectory;
  DiscreteTrajectory<Navigation> descending_trajectory;
  // The so-called North is orthogonal to the plane of the trajectory.
  ComputeNodes(trajectory_in_navigation->Begin(),
               trajectory_in_navigation->End(),
               Vector<double, Navigation>({0, 0, 1}),
               ascending_trajectory,
               descending_trajectory);
  ascending = RenderNavigationTrajectoryInWorld(ascending_trajectory.Begin(),
                                                ascending_trajectory.End(),
                                                sun_world_position);
  descending = RenderNavigationTrajectoryInWorld(descending_trajectory.Begin(),
                                                 descending_trajectory.End(),
                                                 sun_world_position);
}

void Plugin::SetPredictionLength(Time const& t) {
  prediction_length_ = t;
}

void Plugin::SetPredictionAdaptiveStepParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        prediction_adaptive_step_parameters) {
  prediction_parameters_ = prediction_adaptive_step_parameters;
  for (auto const& pair : vessels_) {
    not_null<std::unique_ptr<Vessel>> const& vessel = pair.second;
    vessel->set_prediction_adaptive_step_parameters(prediction_parameters_);
  }
}

bool Plugin::HasVessel(GUID const& vessel_guid) const {
  return Contains(vessels_, vessel_guid);
}

not_null<Vessel*> Plugin::GetVessel(GUID const& vessel_guid) const {
  CHECK(!initializing_);
  return find_vessel_by_guid_or_die(vessel_guid).get();
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

void Plugin::SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) {
  plotting_frame_ = std::move(plotting_frame);
}

not_null<NavigationFrame const*> Plugin::GetPlottingFrame() const {
  return plotting_frame_.get();
}

void Plugin::SetTargetVessel(GUID const& vessel_guid,
                             Index const reference_body_index) {
  not_null<Celestial const*> celestial =
      FindOrDie(celestials_, reference_body_index).get();
  if (!target_ || target_->vessel->guid() != vessel_guid ||
      target_->celestial != celestial) {
    target_.emplace(find_vessel_by_guid_or_die(vessel_guid).get(),
                    ephemeris_.get(),
                    celestial);
  }
  target_->vessel->UpdatePrediction(current_time_ + prediction_length_);
}

void Plugin::ClearTargetVessel() {
  target_ = std::experimental::nullopt;
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
      plugin_->ephemeris_->Prolong(current_time);

      AffineMap<Barycentric, Navigation, Length, OrthogonalMap> const
          barycentric_to_navigation =
              plugin_->plotting_frame_->ToThisFrameAtTime(current_time)
                  .rigid_transformation();
      Position<Navigation> const q_in_navigation =
          (barycentric_to_navigation *
           plugin_->WorldToBarycentric(sun_world_position_))(q);

      OrthogonalMap<RightHandedNavball, Barycentric>
          right_handed_navball_to_barycentric =
              barycentric_right_handed_field_ == nullptr
                  ? plugin_->plotting_frame_->FromThisFrameAtTime(current_time)
                            .orthogonal_map() *
                        navigation_right_handed_field_
                            ->FromThisFrame(q_in_navigation)
                            .Forget()
                  : barycentric_right_handed_field_
                        ->FromThisFrame(
                            plugin_->WorldToBarycentric(sun_world_position_)(q))
                        .Forget();

      // KSP's navball has x west, y up, z south.
      // We want x north, y east, z down.
      OrthogonalMap<Navball, World> const orthogonal_map =
          plugin_->BarycentricToWorld() * right_handed_navball_to_barycentric *
          Permutation<World, RightHandedNavball>(
              Permutation<World, RightHandedNavball>::XZY).Forget() *
          Rotation<Navball, World>(π / 2 * Radian,
                                   Bivector<double, World>({0, 1, 0}),
                                   DefinesFrame<Navball>()).Forget();
      CHECK(orthogonal_map.Determinant().Positive());
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
      dynamic_cast<BodySurfaceDynamicFrame<Barycentric, Navigation>*>(
          &*plotting_frame_);
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
  return FromVesselFrenetFrame(*find_vessel_by_guid_or_die(vessel_guid),
                               Vector<double, Frenet<Navigation>>({1, 0, 0}));
}

Vector<double, World> Plugin::VesselNormal(GUID const& vessel_guid) const {
  return FromVesselFrenetFrame(*find_vessel_by_guid_or_die(vessel_guid),
                               Vector<double, Frenet<Navigation>>({0, 1, 0}));
}

Vector<double, World> Plugin::VesselBinormal(GUID const& vessel_guid) const {
  return FromVesselFrenetFrame(*find_vessel_by_guid_or_die(vessel_guid),
                               Vector<double, Frenet<Navigation>>({0, 0, 1}));
}

Velocity<World> Plugin::VesselVelocity(GUID const& vessel_guid) const {
  Vessel const& vessel = *find_vessel_by_guid_or_die(vessel_guid);
  auto const& last = vessel.psychohistory().last();
  Instant const& time = last.time();
  DegreesOfFreedom<Barycentric> const& barycentric_degrees_of_freedom =
      last.degrees_of_freedom();
  DegreesOfFreedom<Navigation> const plotting_frame_degrees_of_freedom =
      plotting_frame_->ToThisFrameAtTime(time)(barycentric_degrees_of_freedom);
  return Identity<WorldSun, World>()(BarycentricToWorldSun()(
      plotting_frame_->FromThisFrameAtTime(time).orthogonal_map()(
          plotting_frame_degrees_of_freedom.velocity())));
}

AffineMap<Barycentric, World, Length, OrthogonalMap> Plugin::BarycentricToWorld(
    Position<World> const& sun_world_position) const {
  return AffineMap<Barycentric, World, Length, OrthogonalMap>(
      sun_->current_position(current_time_),
      sun_world_position,
      BarycentricToWorld());
}

OrthogonalMap<Barycentric, World> Plugin::BarycentricToWorld() const {
  return OrthogonalMap<WorldSun, World>::Identity() * BarycentricToWorldSun();
}

OrthogonalMap<Barycentric, WorldSun> Plugin::BarycentricToWorldSun() const {
  return sun_looking_glass.Inverse().Forget() * PlanetariumRotation().Forget();
}

AffineMap<World, Barycentric, Length, OrthogonalMap> Plugin::WorldToBarycentric(
    Position<World> const& sun_world_position) const {
  return AffineMap<World, Barycentric, Length, OrthogonalMap>(
      sun_world_position,
      sun_->current_position(current_time_),
      WorldToBarycentric());
}

OrthogonalMap<World, Barycentric> Plugin::WorldToBarycentric() const {
  return BarycentricToWorld().Inverse();
}

Instant Plugin::GameEpoch() const {
  return game_epoch_;
}

Instant Plugin::CurrentTime() const {
  return current_time_;
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
  }
  std::map<not_null<Vessel const*>, GUID const> vessel_to_guid;
  for (auto const& pair : vessels_) {
    std::string const& guid = pair.first;
    not_null<Vessel*> const vessel = pair.second.get();
    vessel_to_guid.emplace(vessel, guid);
    auto* const vessel_message = message->add_vessel();
    vessel_message->set_guid(guid);
    vessel->WriteToMessage(vessel_message->mutable_vessel());
    Index const parent_index = FindOrDie(celestial_to_index, vessel->parent());
    vessel_message->set_parent_index(parent_index);
    vessel_message->set_loaded(Contains(loaded_vessels_, vessel));
    vessel_message->set_kept(Contains(kept_vessels_, vessel));
  }
  for (auto const& pair : part_id_to_vessel_) {
    PartId const part_id = pair.first;
    not_null<Vessel*> const vessel = pair.second;
    (*message->mutable_part_id_to_vessel())[part_id] = vessel_to_guid[vessel];
  }

  ephemeris_->WriteToMessage(message->mutable_ephemeris());

  history_parameters_.WriteToMessage(message->mutable_history_parameters());
  prolongation_parameters_.WriteToMessage(
      message->mutable_prolongation_parameters());
  prediction_parameters_.WriteToMessage(
      message->mutable_prediction_parameters());

  planetarium_rotation_.WriteToMessage(message->mutable_planetarium_rotation());
  game_epoch_.WriteToMessage(message->mutable_game_epoch());
  current_time_.WriteToMessage(message->mutable_current_time());
  Index const sun_index = FindOrDie(celestial_to_index, sun_);
  message->set_sun_index(sun_index);
  plotting_frame_->WriteToMessage(message->mutable_plotting_frame());

  for (auto const& pile_up : pile_ups_) {
    pile_up.WriteToMessage(message->add_pile_up());
  }

  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

not_null<std::unique_ptr<Plugin>> Plugin::ReadFromMessage(
    serialization::Plugin const& message) {
  LOG(INFO) << __FUNCTION__;

  auto const history_parameters =
      Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
          message.history_parameters());
  auto const prolongation_parameters =
      Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
          message.prolongation_parameters());
  auto const prediction_parameters =
      Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
          message.prediction_parameters());
  not_null<std::unique_ptr<Plugin>> plugin =
      std::unique_ptr<Plugin>(new Plugin(history_parameters,
                                         prolongation_parameters,
                                         prediction_parameters));

  plugin->ephemeris_ =
      Ephemeris<Barycentric>::ReadFromMessage(message.ephemeris());
  ReadCelestialsFromMessages(*plugin->ephemeris_,
                             message.celestial(),
                             plugin->celestials_);

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
    auto const inserted =
        plugin->vessels_.emplace(vessel_message.guid(), std::move(vessel));
    CHECK(inserted.second);
  }

  for (auto const& pair : message.part_id_to_vessel()) {
    PartId const part_id = pair.first;
    GUID const guid = pair.second;
    auto const& vessel = FindOrDie(plugin->vessels_, guid);
    plugin->part_id_to_vessel_.emplace(part_id, vessel.get());
  }

  plugin->game_epoch_ = Instant::ReadFromMessage(message.game_epoch());
  plugin->current_time_ = Instant::ReadFromMessage(message.current_time());
  plugin->planetarium_rotation_ =
      Angle::ReadFromMessage(message.planetarium_rotation());

  plugin->sun_ = FindOrDie(plugin->celestials_, message.sun_index()).get();
  plugin->main_body_ = plugin->sun_->body();
  plugin->UpdatePlanetariumRotation();

  std::unique_ptr<NavigationFrame> plotting_frame =
      NavigationFrame::ReadFromMessage(plugin->ephemeris_.get(),
                                       message.plotting_frame());
  plugin->SetPlottingFrame(std::move(plotting_frame));

  // Note that for proper deserialization of parts this list must be
  // reconstructed in its original order.
  for (auto const& pile_up_message : message.pile_up()) {
    plugin->pile_ups_.push_back(PileUp::ReadFromMessage(
        pile_up_message,
        [&part_id_to_vessel = plugin->part_id_to_vessel_](
            PartId const part_id) {
          not_null<Vessel*> const vessel = part_id_to_vessel.at(part_id);
          not_null<Part*> const part = vessel->part(part_id);
          return part;
        },
        plugin->ephemeris_.get()));
  }

  // Now fill the containing pile-up of all the parts.
  for (auto const& vessel_message : message.vessel()) {
    GUID const guid = vessel_message.guid();
    auto const& vessel = FindOrDie(plugin->vessels_, guid);
    vessel->FillContainingPileUpsFromMessage(vessel_message.vessel(),
                                             &plugin->pile_ups_);
  }

  plugin->initializing_.Flop();
  return std::move(plugin);
}

std::unique_ptr<Ephemeris<Barycentric>> Plugin::NewEphemeris(
    std::vector<
        not_null<std::unique_ptr<RotatingBody<Barycentric> const>>>&& bodies,
    std::vector<DegreesOfFreedom<Barycentric>> const& initial_state,
    Instant const& initial_time,
    Length const& fitting_tolerance,
    Ephemeris<Barycentric>::FixedStepParameters const& parameters) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> massive_bodies;
  for (auto& body : bodies) {
    massive_bodies.emplace_back(std::move(body));
  }
  return std::make_unique<Ephemeris<Barycentric>>(std::move(massive_bodies),
                                                  initial_state,
                                                  initial_time,
                                                  fitting_tolerance,
                                                  parameters);
}

Plugin::Plugin(
    Ephemeris<Barycentric>::FixedStepParameters const& history_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        prolongation_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters const& prediction_parameters)
    : history_parameters_(history_parameters),
      prolongation_parameters_(prolongation_parameters),
      prediction_parameters_(prediction_parameters) {}

void Plugin::InitializeEphemerisAndSetCelestialTrajectories() {
  std::vector<
      not_null<std::unique_ptr<RotatingBody<Barycentric> const>>> bodies;
  std::vector<DegreesOfFreedom<Barycentric>> initial_state;
  for (auto& pair : absolute_initialization_->bodies) {
    auto& body = pair.second;
    bodies.emplace_back(std::move(body));
  }
  for (auto const& pair : absolute_initialization_->initial_state) {
    auto const& degrees_of_freedom = pair.second;
    initial_state.emplace_back(degrees_of_freedom);
  }
  absolute_initialization_ = std::experimental::nullopt;
  ephemeris_ = NewEphemeris(std::move(bodies),
                            initial_state,
                            current_time_,
                            default_ephemeris_fitting_tolerance,
                            DefaultEphemerisParameters());
  for (auto const& pair : celestials_) {
    auto& celestial = *pair.second;
    celestial.set_trajectory(ephemeris_->trajectory(celestial.body()));
  }

  // This would use NewBodyCentredNonRotatingNavigationFrame, but we don't have
  // the sun's index at hand.
  // TODO(egg): maybe these functions should take |Celestial*|s, and we should
  // then export |FindOrDie(celestials_, _)|.
  SetPlottingFrame(
      make_not_null_unique<
          BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>>(
              ephemeris_.get(),
              sun_->body()));
}

not_null<std::unique_ptr<Vessel>> const& Plugin::find_vessel_by_guid_or_die(
    GUID const& vessel_guid) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel_guid);
  VLOG_AND_RETURN(1, FindOrDie(vessels_, vessel_guid));
}

// The map between the vector spaces of |Barycentric| and |AliceSun| at
// |current_time_|.
Rotation<Barycentric, AliceSun> const& Plugin::PlanetariumRotation() const {
  return *cached_planetarium_rotation_;
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

Vector<double, World> Plugin::FromVesselFrenetFrame(
    Vessel const& vessel,
    Vector<double, Frenet<Navigation>> const& vector) const {
  auto const& last = vessel.psychohistory().last();
  Instant const& time = last.time();
  DegreesOfFreedom<Barycentric> const& degrees_of_freedom =
      last.degrees_of_freedom();
  auto const from_frenet_frame_to_navigation_frame =
      plotting_frame_->FrenetFrame(
          time,
          plotting_frame_->ToThisFrameAtTime(time)(degrees_of_freedom));

  // The given |vector| in the Frenet frame of the vessel's free-falling
  // trajectory in the given |navigation_frame|, converted to |WorldSun|
  // coordinates.
  return Identity<WorldSun, World>()(
      BarycentricToWorldSun()(
          plotting_frame_->FromThisFrameAtTime(time).orthogonal_map()(
              from_frenet_frame_to_navigation_frame(vector))));
}

template<typename T>
void Plugin::ReadCelestialsFromMessages(
  Ephemeris<Barycentric> const& ephemeris,
  google::protobuf::RepeatedPtrField<T> const& celestial_messages,
  IndexToOwnedCelestial& celestials) {
  auto const& bodies = ephemeris.bodies();
  auto bodies_it = bodies.begin();
  for (auto const& celestial_message : celestial_messages) {
    auto const inserted = celestials.emplace(
        celestial_message.index(),
        make_not_null_unique<Celestial>(
            dynamic_cast_not_null<RotatingBody<Barycentric> const*>(
                *bodies_it)));
    CHECK(inserted.second);
    inserted.first->second->set_trajectory(ephemeris.trajectory(*bodies_it));
    ++bodies_it;
  }
  CHECK_EQ(bodies.end() - bodies.begin(), bodies_it - bodies.begin());
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

std::uint64_t Plugin::FingerprintCelestialJacobiKeplerian(
    Index const celestial_index,
    std::experimental::optional<Index> const& parent_index,
    std::experimental::optional<physics::KeplerianElements<Barycentric>> const&
        keplerian_elements,
    RotatingBody<Barycentric> const& body) {
  serialization::CelestialJacobiKeplerian message;
  message.set_celestial_index(celestial_index);
  if (parent_index) {
    message.set_parent_index(*parent_index);
  }
  if (keplerian_elements) {
    keplerian_elements->WriteToMessage(message.mutable_keplerian_elements());
  }
  body.WriteToMessage(message.mutable_body());

  const std::string serialized = message.SerializeAsString();
  return Fingerprint2011(serialized.c_str(), serialized.size());
}

not_null<std::unique_ptr<DiscreteTrajectory<Navigation>>>
Plugin::RenderBarycentricTrajectoryInNavigation(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end) const {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<Navigation>>();

  NavigationFrame& plotting_frame =
      target_ ? *target_->target_frame : *plotting_frame_;

  if (target_ && !begin.trajectory()->Empty() &&
      (target_->vessel->prediction().Empty() ||
       begin.trajectory()->last().time() >
           target_->vessel->prediction().last().time())) {
    // NOTE(egg): this is an ugly hack to try to get a long enough trajectory
    // while retaining a timeout.
    auto parameters = target_->vessel->prediction_adaptive_step_parameters();
    parameters.set_max_steps(begin.trajectory()->Size());
    target_->vessel->set_prediction_adaptive_step_parameters(parameters);
    target_->vessel->UpdatePrediction(current_time_ + prediction_length_);
  }

  for (auto it = begin; it != end; ++it) {
    if (target_) {
      if (it.time() < target_->vessel->prediction().t_min()) {
        continue;
      } else if (it.time() > target_->vessel->prediction().t_max()) {
        break;
      }
    }
    trajectory->Append(
        it.time(),
        plotting_frame.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
  }
  VLOG(1) << "Returning a " << trajectory->Size() << "-point trajectory";
  return trajectory;
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
Plugin::RenderNavigationTrajectoryInWorld(
    DiscreteTrajectory<Navigation>::Iterator const& begin,
    DiscreteTrajectory<Navigation>::Iterator const& end,
    Position<World> const& sun_world_position) const {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<World>>();

  NavigationFrame& plotting_frame =
      target_ ? *target_->target_frame : *plotting_frame_;

  RigidMotion<Navigation, World> from_navigation_frame_to_world_at_current_time(
      /*rigid_transformation=*/BarycentricToWorld(sun_world_position) *
          plotting_frame.FromThisFrameAtTime(current_time_)
              .rigid_transformation(),
      AngularVelocity<Navigation>{},
      Velocity<Navigation>{});
  for (auto it = begin; it != end; ++it) {
    DegreesOfFreedom<Navigation> const& navigation_degrees_of_freedom =
        it.degrees_of_freedom();
    DegreesOfFreedom<World> const world_degrees_of_freedom =
        from_navigation_frame_to_world_at_current_time(
            navigation_degrees_of_freedom);
    trajectory->Append(it.time(), world_degrees_of_freedom);
  }
  VLOG(1) << "Returning a " << trajectory->Size() << "-point trajectory";
  return trajectory;
}

void Plugin::AddPart(not_null<Vessel*> const vessel,
                     PartId const part_id,
                     std::string const& name,
                     Mass const mass,
                     DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  std::map<PartId, not_null<Vessel*>>::iterator it;
  bool emplaced;
  std::tie(it, emplaced) = part_id_to_vessel_.emplace(part_id, vessel);
  CHECK(emplaced) << NAMED(part_id);
  auto deletion_callback = [it, &map = part_id_to_vessel_] {
    map.erase(it);
  };
  auto part = make_not_null_unique<Part>(part_id,
                                         name,
                                         mass,
                                         degrees_of_freedom,
                                         std::move(deletion_callback));
  vessel->AddPart(std::move(part));
}

bool Plugin::is_loaded(not_null<Vessel*> vessel) const {
  return Contains(loaded_vessels_, vessel);
}

Plugin::Target::Target(not_null<Vessel*> vessel,
                       not_null<Ephemeris<Barycentric> const*> ephemeris,
                       not_null<Celestial const*> const celestial)
    : vessel(vessel),
      celestial(celestial),
      target_frame(
          make_not_null_unique<
              BodyCentredBodyDirectionDynamicFrame<Barycentric, Navigation>>(
              ephemeris,
              [this]() -> auto& { return this->vessel->prediction(); },
              celestial->body())) {}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
