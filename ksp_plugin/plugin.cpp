#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <set>

#include "base/map_util.hpp"
#include "base/not_null.hpp"
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

namespace principia {
namespace ksp_plugin {

using base::FindOrDie;
using base::make_not_null_unique;
using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Identity;
using geometry::Normalize;
using geometry::Permutation;
using geometry::Sign;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::Force;
using si::Radian;

namespace {

// The map between the vector spaces of |World| and |AliceWorld|.
Permutation<World, AliceWorld> const kWorldLookingGlass(
    Permutation<World, AliceWorld>::CoordinatePermutation::XZY);

// The map between the vector spaces of |WorldSun| and |AliceSun|.
Permutation<WorldSun, AliceSun> const kSunLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

}  // namespace

Plugin::Plugin(Instant const& initial_time,
               Index const sun_index,
               GravitationalParameter const& sun_gravitational_parameter,
               Angle const& planetarium_rotation)
    : bubble_(make_not_null_unique<PhysicsBubble>()),
      bodies_(std::make_unique<CelestialToMassiveBody>()),
      initial_state_(std::make_unique<CelestialToDegreesOfFreedom>()),
      history_integrator_(
          McLachlanAtela1992Order5Optimal<Position<Barycentric>>()),
      prolongation_integrator_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>()),
      prediction_integrator_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>()),
      planetarium_rotation_(planetarium_rotation),
      current_time_(initial_time) {
  auto sun_body = std::make_unique<MassiveBody>(sun_gravitational_parameter);
  auto const inserted = celestials_.emplace(sun_index, sun_body.get());
  sun_ = inserted.first->second.get();
  bodies_->emplace(sun_, std::move(sun_body));
  initial_state_->emplace(std::piecewise_construct,
                          std::forward_as_tuple(sun_),
                          std::forward_as_tuple(Position<Barycentric>(),
                                                Velocity<Barycentric>()));
}

void Plugin::InsertCelestial(
    Index const celestial_index,
    GravitationalParameter const& gravitational_parameter,
    Index const parent_index,
    RelativeDegreesOfFreedom<AliceSun> const& from_parent) {
  CHECK(initializing_) << "Celestial bodies should be inserted before the end "
                       << "of initialization";
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  auto const body = std::make_unique<MassiveBody>(gravitational_parameter);
  auto const inserted = celestials_.emplace(celestial_index, body.get());
  CHECK(inserted.second) << "Body already exists at index " << celestial_index;
  LOG(INFO) << "Initial |{orbit.pos, orbit.vel}| for celestial at index "
            << celestial_index << ": " << from_parent;
  auto const relative =
      PlanetariumRotation().Inverse()(from_parent);
  LOG(INFO) << "In barycentric coordinates: " << relative;
  not_null<Celestial*> const celestial = inserted.first->second.get();
  bodies_->emplace(celestial, std::move(body));
  celestial->set_parent(parent);
  DegreesOfFreedom<Barycentric> const& parent_degrees_of_freedom =
      FindOrDie(*initial_state_, parent);
  initial_state_->emplace(celestial,
                          parent_degrees_of_freedom + relative);
}

void Plugin::EndInitialization() {
  initializing_.Flop();
}

void Plugin::UpdateCelestialHierarchy(Index const celestial_index,
                                      Index const parent_index) const {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(celestial_index) << '\n' << NAMED(parent_index);
  CHECK(!initializing_);
  FindOrDie(celestials_, celestial_index)->set_parent(
      FindOrDie(celestials_, parent_index).get());
}

bool Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                Index const parent_index) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(vessel_guid) << '\n' << NAMED(parent_index);
  CHECK(!initializing_);
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  auto inserted = vessels_.emplace(vessel_guid,
                                   make_not_null_unique<Vessel>(parent));
  not_null<Vessel*> const vessel = inserted.first->second.get();
  kept_vessels_.emplace(vessel);
  vessel->set_parent(parent);
  LOG_IF(INFO, inserted.second) << "Inserted vessel with GUID " << vessel_guid
                                << " at " << vessel;
  VLOG(1) << "Parent of vessel with GUID " << vessel_guid <<" is at index "
          << parent_index;
  return inserted.second;
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
  vessel->CreateProlongation(
      current_time_,
      vessel->parent()->trajectory().EvaluateDegreesOfFreedom(
          current_time_,
          vessel->parent()->current_time_hint()) + relative);
  auto const inserted = unsynchronized_vessels_.emplace(vessel.get());
  CHECK(inserted.second);
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(t) << '\n' << NAMED(planetarium_rotation);
  CHECK(!initializing_);
  CHECK_GT(t, current_time_);
  CleanUpVessels();
  bubble_->Prepare(BarycentricToWorldSun(), current_time_, t);
  if (HistoryTime() + Δt_ < t) {
    // The histories are far enough behind that we can advance them at least one
    // step and reset the prolongations.
    EvolveHistories(t);
    // TODO(egg): I think |!bubble_->empty()| => |has_dirty_vessels()|.
    if (has_unsynchronized_vessels() ||
        has_dirty_vessels() ||
        !bubble_->empty()) {
      SynchronizeNewVesselsAndCleanDirtyVessels();
    }
    ResetProlongations();
  }
  EvolveProlongationsAndBubble(t);
  VLOG(1) << "Time has been advanced" << '\n'
          << "from : " << current_time_ << '\n'
          << "to   : " << t;
  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
  UpdatePredictions();
}

void Plugin::ForgetAllHistoriesBefore(Instant const& t) const {
  CHECK(!initializing_);
  CHECK_LT(t, HistoryTime());
  n_body_system_->ForgetBefore(t);
  for (auto const& pair : vessels_) {
    not_null<std::unique_ptr<Vessel>> const& vessel = pair.second;
    // Only forget the synchronized vessels, the others don't have an history.
    if (unsynchronized_vessels_.count(vessel.get()) == 0) {
      vessel->mutable_history()->ForgetBefore(t);
    }
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
      vessel->parent()->prolongation().last().degrees_of_freedom();
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
  Celestial const& celestial = *FindOrDie(celestials_, celestial_index);
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      celestial.prolongation().last().degrees_of_freedom() -
      celestial.parent()->prolongation().last().degrees_of_freedom();
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  VLOG(1) << "Celestial at index " << celestial_index
          << " is at parent degrees of freedom + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

RenderedTrajectory<World> Plugin::RenderedVesselTrajectory(
    GUID const& vessel_guid,
    not_null<RenderingTransforms*> const transforms,
    Position<World> const& sun_world_position) const {
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized());
  VLOG(1) << "Rendering a trajectory for the vessel with GUID " << vessel_guid;
  if (!vessel->is_synchronized()) {
    // TODO(egg): We render neither unsynchronized histories nor prolongations
    // at the moment.
    VLOG(1) << "Returning an empty trajectory";
    return RenderedTrajectory<World>();
  }

  // Compute the apparent trajectory using the given |transforms|.
  return RenderTrajectory(vessel->body(),
                          transforms->first(*vessel,
                                            &MobileInterface::history),
                          transforms,
                          sun_world_position);
}

RenderedTrajectory<World> Plugin::RenderedPrediction(
    not_null<RenderingTransforms*> const transforms,
    Position<World> const& sun_world_position) {
  CHECK(!initializing_);
  if (!HasPredictions()) {
    return RenderedTrajectory<World>();
  }
  RenderedTrajectory<World> result =
      RenderTrajectory(predicted_vessel_->body(),
                       transforms->first_on_or_after(
                           *predicted_vessel_,
                           &MobileInterface::prediction,
                           *predicted_vessel_->prediction().fork_time()),
                       transforms,
                       sun_world_position);
  return result;
}

void Plugin::set_predicted_vessel(GUID const& vessel_guid) {
  clear_predicted_vessel();
  predicted_vessel_ = find_vessel_by_guid_or_die(vessel_guid).get();
}

void Plugin::clear_predicted_vessel() {
  DeletePredictions();
  predicted_vessel_ = nullptr;
}

void Plugin::set_prediction_length(Time const& t) {
  prediction_length_ = t;
}

bool Plugin::has_vessel(GUID const& vessel_guid) const {
  return vessels_.find(vessel_guid) != vessels_.end();
}

not_null<std::unique_ptr<RenderingTransforms>>
Plugin::NewBodyCentredNonRotatingTransforms(
    Index const reference_body_index) const {
  // TODO(egg): this should be const, use a custom comparator in the map.
  not_null<Celestial*> reference_body =
      FindOrDie(celestials_, reference_body_index).get();
  auto transforms = RenderingTransforms::BodyCentredNonRotating(
                        *reference_body,
                        &MobileInterface::prolongation);
  transforms->set_cacheable(&MobileInterface::history);
  return transforms;
}

not_null<std::unique_ptr<RenderingTransforms>>
Plugin::NewBarycentricRotatingTransforms(Index const primary_index,
                                         Index const secondary_index) const {
  // TODO(egg): these should be const, use a custom comparator in the map.
  not_null<Celestial*> primary =
      FindOrDie(celestials_, primary_index).get();
  not_null<Celestial*> secondary =
      FindOrDie(celestials_, secondary_index).get();
  auto transforms = RenderingTransforms::BarycentricRotating(
                        *primary,
                        *secondary,
                        &MobileInterface::prolongation);
  transforms->set_cacheable(&MobileInterface::history);
  return transforms;
}

void Plugin::AddVesselToNextPhysicsBubble(
    GUID const& vessel_guid,
    std::vector<IdAndOwnedPart> parts) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel_guid) << '\n' << NAMED(parts);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  dirty_vessels_.insert(vessel.get());
  bubble_->AddVesselToNext(vessel.get(), std::move(parts));
}

bool Plugin::PhysicsBubbleIsEmpty() const {
  VLOG(1) << __FUNCTION__;
  VLOG_AND_RETURN(1, bubble_->empty());
}

Displacement<World> Plugin::BubbleDisplacementCorrection(
    Position<World> const& sun_world_position) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(sun_world_position);
  VLOG_AND_RETURN(1, bubble_->DisplacementCorrection(BarycentricToWorldSun(),
                                                     *sun_,
                                                     sun_world_position));
}

Velocity<World> Plugin::BubbleVelocityCorrection(
    Index const reference_body_index) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(reference_body_index);
  Celestial const& reference_body =
      *FindOrDie(celestials_, reference_body_index);
  VLOG_AND_RETURN(1, bubble_->VelocityCorrection(BarycentricToWorldSun(),
                                                 reference_body));
}

FrameField<World> Plugin::Navball(
    not_null<RenderingTransforms*> const transforms,
    Position<World> const& sun_world_position) const {
  auto const to_world =
      OrthogonalMap<WorldSun, World>::Identity() * BarycentricToWorldSun();
  auto const positions_from_world =
      AffineMap<World, Barycentric, Length, OrthogonalMap>(
          sun_world_position,
          sun_->prolongation().last().degrees_of_freedom().position(),
          to_world.Inverse());
  return [transforms, to_world, positions_from_world](
      Position<World> const& q) -> Rotation<World, World> {
    // KSP's navball has x west, y up, z south.
    // we want x north, y west, z up.
    auto const orthogonal_map = to_world *
        transforms->coordinate_frame()(positions_from_world(q)).Forget() *
        Permutation<World, Barycentric>(
            Permutation<World, Barycentric>::XZY).Forget() *
        Rotation<World, World>(π / 2 * Radian,
                               Bivector<double, World>({0, 1, 0})).Forget();
    CHECK(orthogonal_map.Determinant().Positive());
    return orthogonal_map.rotation();
  };
}

Vector<double, World> Plugin::VesselTangent(
    GUID const& vessel_guid,
    not_null<RenderingTransforms*> const transforms) const {
  Vessel const& vessel = *find_vessel_by_guid_or_die(vessel_guid);
  auto const actual_it =
      transforms->first_on_or_after(vessel,
                                    &MobileInterface::prolongation,
                                    vessel.prolongation().last().time());
  Trajectory<Rendering> intermediate_trajectory(vessel.body());
  intermediate_trajectory.Append(actual_it.time(),
                                 actual_it.degrees_of_freedom());
  auto const intermediate_it = transforms->second(intermediate_trajectory);
  Trajectory<Barycentric> apparent_trajectory(vessel.body());
  apparent_trajectory.Append(intermediate_it.time(),
                             intermediate_it.degrees_of_freedom());
  return Normalize(
    Identity<WorldSun, World>()(
        BarycentricToWorldSun()(
            apparent_trajectory.last().degrees_of_freedom().velocity())));
}

Instant Plugin::current_time() const {
  return current_time_;
}

void Plugin::WriteToMessage(
    not_null<serialization::Plugin*> const message) const {
  LOG(INFO) << __FUNCTION__;
  CHECK(!initializing_);
  std::map<not_null<Celestial const*>, Index const> celestial_to_index;
  for (auto const& index_celestial : celestials_) {
    celestial_to_index.emplace(index_celestial.second.get(),
                               index_celestial.first);
  }
  for (auto const& index_celestial : celestials_) {
    Index const index = index_celestial.first;
    not_null<Celestial const*> const celestial = index_celestial.second.get();
    auto const celestial_message = message->add_celestial();
    celestial_message->set_index(index);
    celestial->WriteToMessage(celestial_message->mutable_celestial());
    if (celestial->has_parent()) {
      Index const parent_index =
          FindOrDie(celestial_to_index, celestial->parent());
      celestial_message->set_parent_index(parent_index);
    }
  }
  std::map<not_null<Vessel const*>, GUID const> vessel_to_guid;
  for (auto const& guid_vessel : vessels_) {
    std::string const& guid = guid_vessel.first;
    not_null<Vessel*> const vessel = guid_vessel.second.get();
    vessel_to_guid.emplace(vessel, guid);
    auto* const vessel_message = message->add_vessel();
    vessel_message->set_guid(guid);
    vessel->WriteToMessage(vessel_message->mutable_vessel());
    Index const parent_index = FindOrDie(celestial_to_index, vessel->parent());
    vessel_message->set_parent_index(parent_index);
    vessel_message->set_dirty(is_dirty(vessel));
  }

  bubble_->WriteToMessage(
      [&vessel_to_guid](not_null<Vessel const*> const vessel) -> GUID {
        return FindOrDie(vessel_to_guid, vessel);
      },
      message->mutable_bubble());

  planetarium_rotation_.WriteToMessage(message->mutable_planetarium_rotation());
  current_time_.WriteToMessage(message->mutable_current_time());
  Index const sun_index = FindOrDie(celestial_to_index, sun_);
  message->set_sun_index(sun_index);
}

std::unique_ptr<Plugin> Plugin::ReadFromMessage(
    serialization::Plugin const& message) {
  LOG(INFO) << __FUNCTION__;
  IndexToOwnedCelestial celestials;
  for (auto const& celestial_message : message.celestial()) {
    celestials.emplace(
        celestial_message.index(),
        Celestial::ReadFromMessage(celestial_message.celestial()));
  }
  for (auto const& celestial_message : message.celestial()) {
    if (celestial_message.has_parent_index()) {
      not_null<std::unique_ptr<Celestial>> const& celestial =
          FindOrDie(celestials, celestial_message.index());
      not_null<Celestial const*> const parent =
          FindOrDie(celestials, celestial_message.parent_index()).get();
      celestial->set_parent(parent);
    }
  }
  GUIDToOwnedVessel vessels;
  std::set<not_null<Vessel*> const> dirty_vessels;
  for (auto const& vessel_message : message.vessel()) {
    not_null<Celestial const*> const parent =
        FindOrDie(celestials, vessel_message.parent_index()).get();
    not_null<std::unique_ptr<Vessel>> vessel =
        Vessel::ReadFromMessage(vessel_message.vessel(), parent);
    if (vessel_message.dirty()) {
      dirty_vessels.emplace(vessel.get());
    }
    auto const inserted =
        vessels.emplace(vessel_message.guid(), std::move(vessel));
    CHECK(inserted.second);
  }
  not_null<std::unique_ptr<PhysicsBubble>> bubble =
      PhysicsBubble::ReadFromMessage(
          [&vessels](GUID guid) -> not_null<Vessel*> {
            return FindOrDie(vessels, guid).get();
          },
          message.bubble());
  // Can't use |make_unique| here without implementation-dependent friendships.
  return std::unique_ptr<Plugin>(
      new Plugin(std::move(vessels),
                 std::move(celestials),
                 std::move(dirty_vessels),
                 std::move(bubble),
                 Angle::ReadFromMessage(message.planetarium_rotation()),
                 Instant::ReadFromMessage(message.current_time()),
                 message.sun_index()));
}

Plugin::Plugin(GUIDToOwnedVessel vessels,
               IndexToOwnedCelestial celestials,
               std::set<not_null<Vessel*> const> dirty_vessels,
               not_null<std::unique_ptr<PhysicsBubble>> bubble,
               Angle planetarium_rotation,
               Instant current_time,
               Index sun_index)
    : vessels_(std::move(vessels)),
      celestials_(std::move(celestials)),
      dirty_vessels_(std::move(dirty_vessels)),
      bubble_(std::move(bubble)),
      history_integrator_(
          McLachlanAtela1992Order5Optimal<Position<Barycentric>>()),
      prolongation_integrator_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>()),
      prediction_integrator_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>()),
      planetarium_rotation_(planetarium_rotation),
      current_time_(current_time),
      sun_(FindOrDie(celestials_, sun_index).get()) {
  for (auto const& guid_vessel : vessels_) {
    auto const& vessel = guid_vessel.second;
    kept_vessels_.emplace(vessel.get());
    if (!vessel->is_synchronized()) {
      unsynchronized_vessels_.emplace(vessel.get());
    }
  }
  EndInitialization();
}

not_null<std::unique_ptr<Vessel>> const& Plugin::find_vessel_by_guid_or_die(
    GUID const& vessel_guid) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel_guid);
  VLOG_AND_RETURN(1, FindOrDie(vessels_, vessel_guid));
}

bool Plugin::has_dirty_vessels() const {
  return !dirty_vessels_.empty();
}

bool Plugin::has_unsynchronized_vessels() const {
  return !unsynchronized_vessels_.empty();
}

bool Plugin::is_dirty(not_null<Vessel*> const vessel) const {
  return dirty_vessels_.count(vessel) > 0;
}

bool Plugin::has_predicted_vessel() const {
  return predicted_vessel_ != nullptr;
}

bool Plugin::HasPredictions() const {
  if (has_predicted_vessel()) {
    bool const has_prediction = predicted_vessel_->has_prediction();
    for (auto const& pair : celestials_) {
      not_null<std::unique_ptr<Celestial>> const& celestial = pair.second;
      CHECK_EQ(celestial->has_prediction(), has_prediction);
    }
    return has_prediction;
  } else {
    return false;
  }
}

void Plugin::DeletePredictions() {
  if (HasPredictions()) {
    predicted_vessel_->DeletePrediction();
    for (auto const& pair : celestials_) {
      not_null<std::unique_ptr<Celestial>> const& celestial = pair.second;
      celestial->DeletePrediction();
    }
  }
}

Instant const& Plugin::HistoryTime() const {
  return sun_->history().last().time();
}

// The map between the vector spaces of |Barycentric| and |AliceSun| at
// |current_time_|.
Rotation<Barycentric, AliceSun> Plugin::PlanetariumRotation() const {
  return Rotation<Barycentric, AliceSun>(
      planetarium_rotation_,
      Bivector<double, Barycentric>({0, 0, -1}));
}

OrthogonalMap<Barycentric, WorldSun> Plugin::BarycentricToWorldSun() const {
  return kSunLookingGlass.Inverse().Forget() * PlanetariumRotation().Forget();
}

void Plugin::CleanUpVessels() {
  VLOG(1) <<  __FUNCTION__;
  // Remove the vessels which were not updated since last time.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    // While we're going over the vessels, check invariants.
    CheckVesselInvariants(it);
    not_null<Vessel*> const vessel = it->second.get();
    // Now do the cleanup.
    if (kept_vessels_.erase(vessel)) {
      ++it;
    } else {
      LOG(INFO) << "Removing vessel with GUID " << it->first;
      // Since we are going to delete the vessel, we must remove it from
      // |new_vessels| if it's there.
      if (unsynchronized_vessels_.erase(vessel)) {
        LOG(INFO) << "Vessel had not been synchronized";
      }
      if (dirty_vessels_.erase(vessel)) {
        LOG(INFO) << "Vessel was dirty";
      }
      // |std::map::erase| invalidates its parameter so we post-increment.
      vessels_.erase(it++);
    }
  }
}

void Plugin::CheckVesselInvariants(
    GUIDToOwnedVessel::const_iterator const it) const {
  not_null<std::unique_ptr<Vessel>> const& vessel = it->second;
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << it->first
                                  << " was not given an initial state";
  // TODO(egg): At the moment, if a vessel is inserted when
  // |current_time_ == HistoryTime()| (that only happens before the first call
  // to |AdvanceTime|) its first step is unsynchronized. This is convenient to
  // test code paths, but it means the invariant is GE, rather than GT.
  CHECK_GE(vessel->prolongation().last().time(), HistoryTime());
  if (unsynchronized_vessels_.count(vessel.get()) > 0) {
    CHECK(!vessel->is_synchronized());
  } else {
    CHECK(vessel->is_synchronized());
    CHECK_EQ(vessel->history().last().time(), HistoryTime());
  }
}

void Plugin::EvolveHistories(Instant const& t) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(t);
  // Integration with a constant step.
  Ephemeris<Barycentric>::Trajectories trajectories;
  // NOTE(egg): This may be too large, vessels that are not new and in the
  // physics bubble or dirty will not be added.
  trajectories.reserve(vessels_.size() - unsynchronized_vessels_.size());
  for (auto const& pair : vessels_) {
    not_null<Vessel*> const vessel = pair.second.get();
    if (vessel->is_synchronized() &&
        !bubble_->contains(vessel) &&
        !is_dirty(vessel)) {
      trajectories.push_back(vessel->mutable_history());
    }
  }
  VLOG(1) << "Starting the evolution of the histories" << '\n'
          << "from : " << HistoryTime();
  // We integrate until at least |t - Δt_|, and therefore until at most
  // |t|.
  n_body_system_->FlowWithFixedStep(trajectories,
                                    Δt_,
                                    t - Δt_);
  CHECK_GE(HistoryTime(), current_time_);
  VLOG(1) << "Evolved the histories" << '\n'
          << "to   : " << HistoryTime();
}

void Plugin::SynchronizeNewVesselsAndCleanDirtyVessels() {
  VLOG(1) << __FUNCTION__;
  Ephemeris<Barycentric>::Trajectories trajectories;
  trajectories.reserve(unsynchronized_vessels_.size() + bubble_->size());
  for (not_null<Vessel*> const vessel : unsynchronized_vessels_) {
    if (!bubble_->contains(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  for (not_null<Vessel*> const vessel : dirty_vessels_) {
    if (!bubble_->contains(vessel) && vessel->is_synchronized()) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  if (!bubble_->empty()) {
    trajectories.push_back(bubble_->mutable_centre_of_mass_trajectory());
  }
  VLOG(1) << "Starting the synchronization of the new vessels"
          << (bubble_->empty() ? "" : " and of the bubble");
  for (auto const& trajectory : trajectories) {
    n_body_system_->FlowWithAdaptiveStep(trajectory,
                                         prolongation_length_tolerance,
                                         prolongation_speed_tolerance,
                                         prolongation_integrator_,
                                         HistoryTime());
  }
  if (!bubble_->empty()) {
    SynchronizeBubbleHistories();
  }
  for (not_null<Vessel*> const vessel : unsynchronized_vessels_) {
    CHECK(!bubble_->contains(vessel));
    vessel->CreateHistoryAndForkProlongation(
        HistoryTime(),
        vessel->prolongation().last().degrees_of_freedom());
    dirty_vessels_.erase(vessel);
  }
  unsynchronized_vessels_.clear();
  for (not_null<Vessel*> const vessel : dirty_vessels_) {
    CHECK(!bubble_->contains(vessel));
    vessel->mutable_history()->Append(
        HistoryTime(),
        vessel->prolongation().last().degrees_of_freedom());
  }
  dirty_vessels_.clear();
  VLOG(1) << "Synchronized the new vessels"
          << (bubble_->empty() ? "" : " and the bubble");
}

void Plugin::SynchronizeBubbleHistories() {
  VLOG(1) << __FUNCTION__;
  DegreesOfFreedom<Barycentric> const& centre_of_mass =
      bubble_->centre_of_mass_trajectory().last().degrees_of_freedom();
  for (not_null<Vessel*> const vessel : bubble_->vessels()) {
    RelativeDegreesOfFreedom<Barycentric> const& from_centre_of_mass =
        bubble_->from_centre_of_mass(vessel);
    if (vessel->is_synchronized()) {
      vessel->mutable_history()->Append(
          HistoryTime(),
          centre_of_mass + from_centre_of_mass);
    } else {
      vessel->CreateHistoryAndForkProlongation(
          HistoryTime(),
          centre_of_mass + from_centre_of_mass);
      CHECK(unsynchronized_vessels_.erase(vessel));
    }
    CHECK(dirty_vessels_.erase(vessel));
  }
}

void Plugin::ResetProlongations() {
  VLOG(1) << __FUNCTION__;
  for (auto const& pair : vessels_) {
    not_null<std::unique_ptr<Vessel>> const& vessel = pair.second;
    vessel->ResetProlongation(HistoryTime());
  }
  for (auto const& pair : celestials_) {
    not_null<std::unique_ptr<Celestial>> const& celestial = pair.second;
    celestial->ResetProlongation(HistoryTime());
  }
  VLOG(1) << "Prolongations have been reset";
}

void Plugin::EvolveProlongationsAndBubble(Instant const& t) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(t);
  Ephemeris<Barycentric>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() -
                       bubble_->number_of_vessels() + bubble_->size());
  for (auto const& pair : vessels_) {
    not_null<Vessel*> const vessel = pair.second.get();
    if (!bubble_->contains(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  if (!bubble_->empty()) {
    trajectories.push_back(bubble_->mutable_centre_of_mass_trajectory());
  }
  VLOG(1) << "Evolving prolongations"
          << (bubble_->empty() ? "" : " and bubble") << '\n'
          << "from : " << trajectories.front()->last().time() << '\n'
          << "to   : " << t;
  for (auto const& trajectory : trajectories) {
    n_body_system_->FlowWithAdaptiveStep(trajectory,
                                         prolongation_length_tolerance,
                                         prolongation_speed_tolerance,
                                         prolongation_integrator_,
                                         t);
  }
  if (!bubble_->empty()) {
    DegreesOfFreedom<Barycentric> const& centre_of_mass =
        bubble_->centre_of_mass_trajectory().last().degrees_of_freedom();
    for (not_null<Vessel*> vessel : bubble_->vessels()) {
      RelativeDegreesOfFreedom<Barycentric> const& from_centre_of_mass =
          bubble_->from_centre_of_mass(vessel);
      vessel->mutable_prolongation()->Append(
          t,
          centre_of_mass + from_centre_of_mass);
    }
  }
}

void Plugin::UpdatePredictions() {
  DeletePredictions();
  if (has_predicted_vessel()) {
    predicted_vessel_->ForkPrediction();
    n_body_system_->FlowWithAdaptiveStep(
        predicted_vessel_->mutable_prediction(),
        prediction_length_tolerance,
        prediction_speed_tolerance,
        prediction_integrator_,
        current_time_ + prediction_length_);
  }
}

RenderedTrajectory<World> Plugin::RenderTrajectory(
    not_null<Body const*> const body,
    Trajectory<Barycentric>::TransformingIterator<Rendering> const& actual_it,
    not_null<RenderingTransforms*> const transforms,
    Position<World> const& sun_world_position) const {
  RenderedTrajectory<World> result;
  auto const to_world =
      AffineMap<Barycentric, World, Length, OrthogonalMap>(
          sun_->prolongation().last().degrees_of_freedom().position(),
          sun_world_position,
          OrthogonalMap<WorldSun, World>::Identity() * BarycentricToWorldSun());

  // First build the trajectory resulting from the first transform.
  Trajectory<Rendering> intermediate_trajectory(body);
  for (auto it = actual_it; !it.at_end(); ++it) {
    intermediate_trajectory.Append(it.time(), it.degrees_of_freedom());
  }

  // Then build the apparent trajectory using the second transform.
  Trajectory<Barycentric> apparent_trajectory(body);
  for (auto intermediate_it = transforms->second(intermediate_trajectory);
       !intermediate_it.at_end();
       ++intermediate_it) {
    apparent_trajectory.Append(intermediate_it.time(),
                               intermediate_it.degrees_of_freedom());
  }

  // Finally use the apparent trajectory to build the result.
  auto initial_it = apparent_trajectory.first();
  if (!initial_it.at_end()) {
    for (auto final_it = initial_it;
         ++final_it, !final_it.at_end();
         initial_it = final_it) {
      result.emplace_back(to_world(initial_it.degrees_of_freedom().position()),
                          to_world(final_it.degrees_of_freedom().position()));
    }
  }
  VLOG(1) << "Returning a " << result.size() << "-segment trajectory";
  return result;
}

}  // namespace ksp_plugin
}  // namespace principia
