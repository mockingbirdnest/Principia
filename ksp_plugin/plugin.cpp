
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <set>

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
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "ksp_plugin/part_subsets.hpp"
#include "physics/barycentric_rotating_dynamic_frame_body.hpp"
#include "physics/body_centred_body_direction_dynamic_frame.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/frame_field.hpp"
#include "physics/rotating_body.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using base::dynamic_cast_not_null;
using base::Error;
using base::FindOrDie;
using base::Fingerprint2011;
using base::FingerprintCat2011;
using base::make_not_null_unique;
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
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using physics::BarycentricRotatingDynamicFrame;
using physics::BodyCentredBodyDirectionDynamicFrame;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::BodySurfaceDynamicFrame;
using physics::CoordinateFrameField;
using physics::DynamicFrame;
using physics::Frenet;
using physics::KeplerianElements;
using physics::RotatingBody;
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

std::uint64_t const ksp_stock_system_fingerprint = 0xD15286A27180CD31u;
std::uint64_t const ksp_fixed_system_fingerprint = 0x648C354716008328u;

// The map between the vector spaces of |WorldSun| and |AliceSun|.
Permutation<WorldSun, AliceSun> const sun_looking_glass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

Ephemeris<Barycentric>::FixedStepParameters DefaultEphemerisParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
             McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
             /*step=*/45 * Minute);
}

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
    not_null<std::unique_ptr<MassiveBody const>> body) {
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

void Plugin::InsertCelestialJacobiKeplerian(
    Index const celestial_index,
    std::experimental::optional<Index> const& parent_index,
    std::experimental::optional<KeplerianElements<Barycentric>> const&
        keplerian_elements,
    not_null<std::unique_ptr<MassiveBody>> body) {
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
  MassiveBody* const unowned_body = body.get();
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
    std::map<not_null<MassiveBody const*>, Index> bodies_to_indices;
    for (auto const& index_body :
         hierarchical_initialization_->indices_to_bodies) {
      bodies_to_indices[index_body.second] = index_body.first;
    }
    auto const parents = std::move(hierarchical_initialization_->parents);
    hierarchical_initialization_ = std::experimental::nullopt;
#if LOG_KSP_SYSTEM
    std::ofstream file;
    if (system_fingerprint == ksp_stock_system_fingerprint) {
      file.open("ksp_stock_system.proto.hex");
    } else if (system_fingerprint == ksp_fixed_system_fingerprint) {
      file.open("ksp_fixed_system.proto.hex");
    } else {
      file.open("unknown_system.proto.hex");
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
      Index const celestial_index = bodies_to_indices[system.bodies[i].get()];
      InsertCelestialAbsoluteCartesian(
          celestial_index,
          FindOrDie(parents, celestial_index),
          system.degrees_of_freedom[i],
          std::move(system.bodies[i]));
    }
#if LOG_KSP_SYSTEM
    file.close();
#endif
  }
  CHECK(absolute_initialization_);
  CHECK_NOTNULL(sun_);
  main_body_ = CHECK_NOTNULL(
      dynamic_cast_not_null<RotatingBody<Barycentric> const*>(sun_->body()));
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
  main_body_ = dynamic_cast_not_null<RotatingBody<Barycentric> const*>(
      FindOrDie(celestials_, index)->body());
  LOG_IF(FATAL, main_body_ == nullptr) << index;
}

Rotation<BodyWorld, World> Plugin::CelestialRotation(
    Index const index) const {
  // |BodyWorld| with its y and z axes swapped (so that z is the polar axis).
  // The basis is right-handed.
  struct BodyFixed;
  Permutation<BodyWorld, BodyFixed> const body_mirror(
      Permutation<BodyWorld, BodyFixed>::XZY);

  auto const& body = dynamic_cast<RotatingBody<Barycentric> const&>(
      *FindOrDie(celestials_, index)->body());

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
  auto const& body = dynamic_cast<RotatingBody<Barycentric> const&>(
      *FindOrDie(celestials_, celestial_index)->body());
  return body.AngleAt(game_epoch_);
}

Time Plugin::CelestialRotationPeriod(Index const celestial_index) const {
  auto const& body = dynamic_cast<RotatingBody<Barycentric> const&>(
      *FindOrDie(celestials_, celestial_index)->body());
  // The result will be negative if the pole is the negative pole
  // (e.g. for Venus).  This is the convention KSP uses for retrograde rotation.
  return 2 * π * Radian / body.angular_frequency();
}

void Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                Index parent_index,
                                bool loaded,
                                bool& inserted) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(vessel_guid) << '\n'
          << NAMED(parent_index) << '\n'
          << NAMED(loaded) << '\n';
  CHECK(!initializing_);
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  auto pair =
      vessels_.emplace(vessel_guid,
                       make_not_null_unique<Vessel>(parent,
                                                    ephemeris_.get(),
                                                    history_parameters_,
                                                    prolongation_parameters_,
                                                    prediction_parameters_));
  auto const& emplacement = pair.first;
  inserted = pair.second;
  not_null<Vessel*> const vessel = emplacement->second.get();
  kept_vessels_.emplace(vessel);
  vessel->set_parent(parent);
  Subset<Vessel>::MakeSingleton(*vessel, vessel);
  LOG_IF(INFO, inserted) << "Inserted vessel "
                            "with GUID "
                         << vessel_guid << " at " << vessel;
  VLOG(1) << "Parent of vessel "
             "with GUID "
          << vessel_guid << " is at index " << parent_index;
}

void Plugin::InsertOrKeepLoadedPart(PartId const part_id,
                                    Mass const& mass,
                                    not_null<Vessel*> const vessel,
                                    bool& inserted) {
  auto it = part_id_to_vessel_.find(part_id);
  inserted = it == part_id_to_vessel_.end();
  if (inserted) {
    auto const pair = part_id_to_vessel_.emplace(part_id, vessel);
    auto const& it = pair.first;
    bool const emplaced = pair.second;
    CHECK(emplaced);
    auto deletion_callback = [it, &map = part_id_to_vessel_] {
      map.erase(it);
    });
    auto part =
        make_not_null_unique<Part>(part_id, mass, std::move(deletion_callback));
    vessel->AddPart(std::move(part));
  } else {
    not_null<Vessel*> current_vessel = it->second;
    if (vessel == current_vessel) {
     vessel->KeepPart(part_id);
    } else {
      FindOrDie(part_id_to_vessel_, part_id) = vessel;
      vessel->AddPart(current_vessel->ExtractPart(part_id));
    }
  }
}

void Plugin::SetVesselStateOffset(
    GUID const& vessel_guid,
    RelativeDegreesOfFreedom<AliceSun> const& from_parent) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(vessel_guid) << '\n' << NAMED(from_parent);
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(!vessel->is_initialized())
      << "Vessel with GUID " << vessel_guid << " already has a trajectory";
  LOG(INFO) << "Initial |{orbit.pos, orbit.vel}| for vessel with GUID "
            << vessel_guid << ": " << from_parent;
  RelativeDegreesOfFreedom<Barycentric> const relative =
      PlanetariumRotation().Inverse()(from_parent);
  LOG(INFO) << "In barycentric coordinates: " << relative;
  ephemeris_->Prolong(current_time_);
  vessel->CreateHistoryAndForkProlongation(
      current_time_,
      vessel->parent()->current_degrees_of_freedom(current_time_) + relative);
}

void Plugin::FreeVesselsAndCollectPileUps() {
  CHECK(!initializing_);
  FreeVessels();
  for (auto const& pair : vessels_) {
    Vessel& vessel = *pair.second;
    Subset<Vessel>::Find(vessel).mutable_properties().Collect(&pile_ups_);
  }
}

void Plugin::AddPileUpToBubble(std::list<PileUp>::iterator const pile_up) {
  pile_ups_in_bubble_.push_back(pile_up);
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(t) << '\n' << NAMED(planetarium_rotation);
  CHECK(!initializing_);
  CHECK_GT(t, current_time_);
  ephemeris_->Prolong(t);

  if (pile_ups_in_bubble_.empty()) {
    bubble_barycentre_ = std::experimental::nullopt;
  } else {
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> bubble_barycentre;
    for (auto const it : pile_ups_in_bubble_) {
      for (not_null<Vessel*> const vessel : it->vessels()) {
        bubble_barycentre.Add(
            vessel->psychohistory().last().degrees_of_freedom(),
            vessel->mass());
      }
    }
    bubble_barycentre_ = bubble_barycentre.Get();
  }
  pile_ups_in_bubble_.clear();

  VLOG(1) << "Time has been advanced" << '\n'
          << "from : " << current_time_ << '\n'
          << "to   : " << t;
  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
}

DegreesOfFreedom<Barycentric> Plugin::GetBubbleBarycentre() const {
  CHECK(bubble_barycentre_);
  return *bubble_barycentre_;
}

DegreesOfFreedom<World> Plugin::CelestialWorldDegreesOfFreedom(
    Index const index) const {
  RigidMotion<Barycentric, World> barycentric_to_world{
      RigidTransformation<Barycentric, World>{GetBubbleBarycentre().position(),
                                              World::origin,
                                              BarycentricToWorld()},
      AngularVelocity<Barycentric>{},
      GetBubbleBarycentre().velocity()};
  return barycentric_to_world(
      FindOrDie(celestials_, index)->trajectory().EvaluateDegreesOfFreedom(
          current_time_, /*hint=*/nullptr));
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
    GUID const& vessel_guid) const {
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << vessel_guid
                                  << " was not given an initial state";
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      vessel->prolongation().last().degrees_of_freedom() -
      vessel->parent()->current_degrees_of_freedom(current_time_);
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  VLOG(1) << "Vessel with GUID " << vessel_guid
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
Plugin::RenderedVesselTrajectory(
    GUID const& vessel_guid,
    Position<World> const& sun_world_position) const {
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized());
  VLOG(1) << "Rendering a trajectory for the vessel with GUID " << vessel_guid;
  return RenderedTrajectoryFromIterators(vessel->history().Begin(),
                                         vessel->history().End(),
                                         sun_world_position);
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>> Plugin::RenderedPrediction(
    GUID const& vessel_guid,
    Position<World> const& sun_world_position) const {
  CHECK(!initializing_);
  Vessel const& vessel = *find_vessel_by_guid_or_die(vessel_guid);
  return RenderedTrajectoryFromIterators(vessel.prediction().Fork(),
                                         vessel.prediction().End(),
                                         sun_world_position);
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
Plugin::RenderedTrajectoryFromIterators(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position) const {
  auto result = make_not_null_unique<DiscreteTrajectory<World>>();

  // Compute the trajectory in the navigation frame.
  DiscreteTrajectory<Navigation> intermediate_trajectory;
  for (auto it = begin; it != end; ++it) {
    intermediate_trajectory.Append(
        it.time(),
        plotting_frame_->ToThisFrameAtTime(it.time())(
            it.degrees_of_freedom()));
  }

  // Render the trajectory at current time in |World|.
  DiscreteTrajectory<Navigation>::Iterator const intermediate_end =
      intermediate_trajectory.End();
  auto from_navigation_frame_to_world_at_current_time =
      BarycentricToWorld(sun_world_position) *
          plotting_frame_->
              FromThisFrameAtTime(current_time_).rigid_transformation();
  for (auto intermediate_it = intermediate_trajectory.Begin();
       intermediate_it != intermediate_end;
       ++intermediate_it) {
    DegreesOfFreedom<Navigation> const navigation_degrees_of_freedom =
        intermediate_it.degrees_of_freedom();
    DegreesOfFreedom<World> const world_degrees_of_freedom =
        DegreesOfFreedom<World>(
            from_navigation_frame_to_world_at_current_time(
                navigation_degrees_of_freedom.position()),
            from_navigation_frame_to_world_at_current_time.linear_map()(
                navigation_degrees_of_freedom.velocity()));
    result->Append(intermediate_it.time(), world_degrees_of_freedom);
  }
  VLOG(1) << "Returning a " << result->Size() << "-point trajectory";
  return result;
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
  ephemeris_->ComputeApsides(FindOrDie(celestials_, celestial_index)->body(),
                             begin, end,
                             apoapsides_trajectory,
                             periapsides_trajectory);
  apoapsides = RenderedTrajectoryFromIterators(apoapsides_trajectory.Begin(),
                                               apoapsides_trajectory.End(),
                                               sun_world_position);
  periapsides = RenderedTrajectoryFromIterators(periapsides_trajectory.Begin(),
                                                periapsides_trajectory.End(),
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
  return vessels_.find(vessel_guid) != vessels_.end();
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
      dynamic_cast_not_null<RotatingBody<Barycentric> const*>(
          reference_body.body()));
}

void Plugin::SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) {
  plotting_frame_ = std::move(plotting_frame);
}

not_null<NavigationFrame const*> Plugin::GetPlottingFrame() const {
  return plotting_frame_.get();
}

void Plugin::ReportCollision(GUID const& vessel1, GUID const& vessel2) const {
  Vessel& v1 = *FindOrDie(vessels_, vessel1);
  Vessel& v2 = *FindOrDie(vessels_, vessel2);
  Subset<Vessel>::Unite(Subset<Vessel>::Find(v1), Subset<Vessel>::Find(v2));
}

std::unique_ptr<FrameField<World, Navball>> Plugin::NavballFrameField(
    Position<World> const& sun_world_position) const {

  struct RightHandedNavball;

  class NavballFrameField : public FrameField<World, Navball> {
   public:
    NavballFrameField(
        not_null<Plugin const*> const plugin,
        base::not_null<
            std::unique_ptr<FrameField<Navigation, RightHandedNavball>>>
            right_handed_navball_field,
        Position<World> const& sun_world_position)
        : plugin_(plugin),
          right_handed_navball_field_(std::move(right_handed_navball_field)),
          sun_world_position_(sun_world_position) {}

    Rotation<Navball, World> FromThisFrame(
        Position<World> const& q) const override {
      Instant const& current_time = plugin_->current_time_;
      plugin_->ephemeris_->Prolong(current_time);

      OrthogonalMap<Navigation, World> const navigation_to_world =
          plugin_->BarycentricToWorld() *
          plugin_->plotting_frame_->FromThisFrameAtTime(current_time)
              .orthogonal_map();

      AffineMap<Barycentric, Navigation, Length, OrthogonalMap> const
          barycentric_to_navigation =
              plugin_->plotting_frame_->ToThisFrameAtTime(current_time)
                  .rigid_transformation();
      Position<Navigation> const q_in_navigation =
          (barycentric_to_navigation *
           plugin_->WorldToBarycentric(sun_world_position_))(q);

      // KSP's navball has x west, y up, z south.
      // We want x north, y east, z down.
      OrthogonalMap<Navball, World> const orthogonal_map =
          navigation_to_world *
          right_handed_navball_field_->FromThisFrame(q_in_navigation).Forget() *
          Permutation<World, RightHandedNavball>(
              Permutation<World, RightHandedNavball>::XZY)
              .Forget() *
          Rotation<Navball, World>(π / 2 * Radian,
                                   Bivector<double, World>({0, 1, 0}),
                                   DefinesFrame<Navball>())
              .Forget();
      CHECK(orthogonal_map.Determinant().Positive());
      return orthogonal_map.rotation();
    }

   private:
    not_null<Plugin const*> const plugin_;
    base::not_null<
        std::unique_ptr<FrameField<Navigation, RightHandedNavball>>> const
        right_handed_navball_field_;
    Position<World> const sun_world_position_;
  };

  return std::make_unique<NavballFrameField>(
             this,
             base::make_not_null_unique<
                 CoordinateFrameField<Navigation, RightHandedNavball>>(),
             sun_world_position);
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

void Plugin::UpdateAllVesselsInPileUps() {
  for (auto& pile_up : pile_ups_) {
    pile_up.UpdateVesselsInPileUpIfUpdated();
  }
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

bool Plugin::MustRotateBodies() const {
  return !is_pre_cardano_;
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
    vessel_message->set_dirty(vessel->is_dirty());
  }

  ephemeris_->WriteToMessage(message->mutable_ephemeris());

  history_parameters_.WriteToMessage(message->mutable_history_parameters());
  prolongation_parameters_.WriteToMessage(
      message->mutable_prolongation_parameters());
  prediction_parameters_.WriteToMessage(
      message->mutable_prediction_parameters());

  planetarium_rotation_.WriteToMessage(message->mutable_planetarium_rotation());
  if (!is_pre_cardano_) {
    // A pre-Cardano save stays pre-Cardano; we cannot pull rotational
    // properties out of thin air.
    game_epoch_.WriteToMessage(message->mutable_game_epoch());
  }
  current_time_.WriteToMessage(message->mutable_current_time());
  Index const sun_index = FindOrDie(celestial_to_index, sun_);
  message->set_sun_index(sun_index);
  plotting_frame_->WriteToMessage(message->mutable_plotting_frame());
  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

not_null<std::unique_ptr<Plugin>> Plugin::ReadFromMessage(
    serialization::Plugin const& message) {
  LOG(INFO) << __FUNCTION__;
  bool const is_pre_bourbaki = message.pre_bourbaki_celestial_size() > 0;
  std::unique_ptr<Ephemeris<Barycentric>> ephemeris;
  IndexToOwnedCelestial celestials;

  if (is_pre_bourbaki) {
    ephemeris = Ephemeris<Barycentric>::ReadFromPreBourbakiMessages(
        message.pre_bourbaki_celestial(),
        fitting_tolerance,
        DefaultEphemerisParameters());
    ReadCelestialsFromMessages(*ephemeris,
                               message.pre_bourbaki_celestial(),
                               celestials);
  } else {
    ephemeris = Ephemeris<Barycentric>::ReadFromMessage(message.ephemeris());
    ReadCelestialsFromMessages(*ephemeris, message.celestial(), celestials);
  }

  GUIDToOwnedVessel vessels;
  for (auto const& vessel_message : message.vessel()) {
    not_null<Celestial const*> const parent =
        FindOrDie(celestials, vessel_message.parent_index()).get();
    not_null<std::unique_ptr<Vessel>> vessel = Vessel::ReadFromMessage(
                                                   vessel_message.vessel(),
                                                   ephemeris.get(),
                                                   parent);
    if (vessel_message.dirty()) {
      vessel->set_dirty();
    }
    Subset<Vessel>::MakeSingleton(*vessel, vessel.get());
    auto const inserted =
        vessels.emplace(vessel_message.guid(), std::move(vessel));
    CHECK(inserted.second);
  }

  Instant const current_time = Instant::ReadFromMessage(message.current_time());

  bool const is_pre_буняковский = !(message.has_history_parameters() &&
                                    message.has_prolongation_parameters() &&
                                    message.has_prediction_parameters());
  auto const history_parameters =
    is_pre_буняковский
        ? DefaultHistoryParameters()
        : Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
              message.history_parameters());
  auto const prolongation_parameters =
    is_pre_буняковский
        ? DefaultProlongationParameters()
        : Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
              message.prolongation_parameters());
  auto const prediction_parameters =
    is_pre_буняковский
        ? DefaultPredictionParameters()
        : Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
              message.prediction_parameters());

  bool const is_pre_cardano = !message.has_game_epoch();
  Instant const game_epoch =
      is_pre_cardano ? astronomy::J2000
                     : Instant::ReadFromMessage(message.game_epoch());

  // Can't use |make_unique| here without implementation-dependent friendships.
  auto plugin = std::unique_ptr<Plugin>(
      new Plugin(std::move(vessels),
                 std::move(celestials),
                 std::move(ephemeris),
                 history_parameters,
                 prolongation_parameters,
                 prediction_parameters,
                 Angle::ReadFromMessage(message.planetarium_rotation()),
                 game_epoch,
                 current_time,
                 message.sun_index(),
                 is_pre_cardano));
  std::unique_ptr<NavigationFrame> plotting_frame =
      NavigationFrame::ReadFromMessage(plugin->ephemeris_.get(),
                                       message.plotting_frame());
  if (plotting_frame == nullptr) {
    // In the pre-Brouwer compatibility case you get a plotting frame centred on
    // the Sun.
    plugin->SetPlottingFrame(
        plugin->NewBodyCentredNonRotatingNavigationFrame(message.sun_index()));
  } else {
    plugin->SetPlottingFrame(std::move(plotting_frame));
  }
  return std::move(plugin);
}

std::unique_ptr<Ephemeris<Barycentric>> Plugin::NewEphemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody const>>>&& bodies,
    std::vector<DegreesOfFreedom<Barycentric>> const& initial_state,
    Instant const& initial_time,
    Length const& fitting_tolerance,
    Ephemeris<Barycentric>::FixedStepParameters const& parameters) {
  return std::make_unique<Ephemeris<Barycentric>>(std::move(bodies),
                                                  initial_state,
                                                  initial_time,
                                                  fitting_tolerance,
                                                  parameters);
}

Plugin::Plugin(
    GUIDToOwnedVessel vessels,
    IndexToOwnedCelestial celestials,
    std::unique_ptr<Ephemeris<Barycentric>> ephemeris,
    Ephemeris<Barycentric>::FixedStepParameters const& history_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        prolongation_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters const& prediction_parameters,
    Angle const& planetarium_rotation,
    Instant const& game_epoch,
    Instant const& current_time,
    Index const sun_index,
    bool const is_pre_cardano)
    : vessels_(std::move(vessels)),
      celestials_(std::move(celestials)),
      ephemeris_(std::move(ephemeris)),
      history_parameters_(history_parameters),
      prolongation_parameters_(prolongation_parameters),
      prediction_parameters_(prediction_parameters),
      planetarium_rotation_(planetarium_rotation),
      game_epoch_(game_epoch),
      current_time_(current_time),
      sun_(FindOrDie(celestials_, sun_index).get()),
      is_pre_cardano_(is_pre_cardano) {
  for (auto const& pair : vessels_) {
    auto const& vessel = pair.second;
    kept_vessels_.emplace(vessel.get());
  }
  if (!is_pre_cardano_) {
    main_body_ = CHECK_NOTNULL(
        dynamic_cast_not_null<RotatingBody<Barycentric> const*>(sun_->body()));
  }
  initializing_.Flop();
}


void Plugin::InitializeEphemerisAndSetCelestialTrajectories() {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
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
                            fitting_tolerance,
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
Rotation<Barycentric, AliceSun> Plugin::PlanetariumRotation() const {
  // The z axis of |PlanetariumFrame| is the pole of |main_body_|, and its x
  // axis is the origin of body rotation (the intersection between the
  // |Barycentric| xy plane and the plane of |main_body_|'s equator, or the y
  // axis of |Barycentric| if they coincide).
  // This can be expressed using Euler angles, see figures 1 and 2 of
  // http://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2009reprint.pdf.
  struct PlanetariumFrame;

  if (is_pre_cardano_) {
    return Rotation<Barycentric, AliceSun>(
        planetarium_rotation_,
        Bivector<double, Barycentric>({0, 0, 1}),
        DefinesFrame<AliceSun>{});
  } else {
    CHECK_NOTNULL(main_body_);
    Rotation<Barycentric, PlanetariumFrame> const to_planetarium(
        π / 2 * Radian + main_body_->right_ascension_of_pole(),
        π / 2 * Radian - main_body_->declination_of_pole(),
        0 * Radian,
        EulerAngles::ZXZ,
        DefinesFrame<PlanetariumFrame>{});
    return Rotation<PlanetariumFrame, AliceSun>(
               planetarium_rotation_,
               Bivector<double, PlanetariumFrame>({0, 0, 1}),
               DefinesFrame<AliceSun>{}) * to_planetarium;
  }
}

void Plugin::FreeVessels() {
  VLOG(1) <<  __FUNCTION__;
  // Remove the vessels which were not updated since last time.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    // While we're going over the vessels, check invariants.
    not_null<Vessel*> const vessel = it->second.get();
    // Now do the cleanup.
    if (kept_vessels_.erase(vessel)) {
      ++it;
    } else {
      LOG(INFO) << "Removing vessel with GUID " << it->first;
      vessel->clear_pile_up();
      it = vessels_.erase(it);
    }
  }
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
        make_not_null_unique<Celestial>(*bodies_it));
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
    MassiveBody const& body) {
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

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
