#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "geometry/identity.hpp"
#include "geometry/permutation.hpp"
#include "geometry/affine_map.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::AffineMap;
using geometry::Bivector;
using geometry::Identity;
using geometry::Permutation;
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

// The map between the vector spaces of |Barycentric| and |WorldSun| at
// |current_time_|.
Rotation<Barycentric, WorldSun> Plugin::PlanetariumRotation() const {
  return Rotation<Barycentric, WorldSun>(
      planetarium_rotation_,
      Bivector<double, Barycentric>({0, 1, 0}));
}

void Plugin::CheckVesselInvariants(
    GUIDToOwnedVessel::const_iterator const it) const {
  Vessel* const vessel = it->second.get();
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << it->first
                                  << " was not given an initial state";
  // TODO(egg): At the moment, if a vessel is inserted when
  // |current_time_ == HistoryTime()| (that only happens before the first call
  // to |AdvanceTime|) its first step is unsynchronized. This is convenient to
  // test code paths, but it means the invariant is GE, rather than GT.
  CHECK_GE(vessel->prolongation().last().time(), HistoryTime());
  if (new_vessels_.count(vessel) > 0) {
    CHECK(!vessel->is_synchronized());
  } else {
    CHECK(vessel->is_synchronized());
    CHECK_EQ(vessel->history().last().time(), HistoryTime());
  }
}

void Plugin::CleanUpVessels() {
  VLOG(1) << "Vessel cleanup";
  // Remove the vessels which were not updated since last time.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    // While we're going over the vessels, check invariants.
    CheckVesselInvariants(it);
    Vessel* const vessel = it->second.get();
    // Now do the cleanup.
    if (kept_.erase(vessel)) {
      ++it;
    } else {
      LOG(INFO) << "Removing vessel with GUID " << it->first;
      // Since we are going to delete the vessel, we must remove it from
      // |new_vessels| if it's there.
      if (new_vessels_.erase(vessel)) {
        LOG(INFO) << "Vessel had not been synchronized";
      }
      // |std::map::erase| invalidates its parameter so we post-increment.
      vessels_.erase(it++);
    }
  }
}

void Plugin::EvolveSynchronizedHistories(Instant const& t) {
  VLOG(1) << "Starting the evolution of the old histories" << '\n'
          << "from : " << HistoryTime();
  // Integration with a constant step.
  NBodySystem<Barycentric>::Trajectories trajectories;
  // NOTE(egg): This may be too large, vessels that are not new and in the
  // physics bubble will not be added.
  trajectories.reserve(vessels_.size() - new_vessels_.size() +
                       celestials_.size());
  for (auto const& pair : celestials_) {
    Celestial* const celestial = pair.second.get();
    trajectories.push_back(celestial->mutable_history());
  }
  for (auto const& pair : vessels_) {
    Vessel* const vessel = pair.second.get();
    if (vessel->is_synchronized() && !IsInPhysicsBubble(vessel)) {
      trajectories.push_back(vessel->mutable_history());
    }
  }
  n_body_system_->Integrate(history_integrator_,  // integrator
                            t,                    // tmax
                            Δt_,                  // Δt
                            0,                    // sampling_period
                            false,                // tmax_is_exact
                            trajectories);        // trajectories
  CHECK_GE(HistoryTime(), current_time_);
  VLOG(1) << "Evolved the old histories" << '\n'
          << "to   : " << HistoryTime();
}

void Plugin::SynchronizeNewHistoriesAndBubble() {
  VLOG(1) << "Starting the synchronization of the new histories";
  NBodySystem<Barycentric>::Trajectories trajectories;
  trajectories.reserve(celestials_.size() + new_vessels_.size() +
                       HavePhysicsBubble() ? 1 : 0);
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (Vessel* const vessel : new_vessels_) {
    if (!IsInPhysicsBubble(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  if (HavePhysicsBubble()) {
    trajectories.push_back(
        current_physics_bubble_->centre_of_mass_trajectory.get());
  }
  n_body_system_->Integrate(prolongation_integrator_,  // integrator
                            HistoryTime(),             // tmax
                            Δt_,                       // Δt
                            0,                         // sampling_period
                            true,                      // tmax_is_exact
                            trajectories);             // trajectories
  if (HavePhysicsBubble()) {
    SynchronizeBubbleHistories();
  }
  for (Vessel* const vessel : new_vessels_) {
    CHECK(!IsInPhysicsBubble(vessel));
    vessel->CreateHistoryAndForkProlongation(
        HistoryTime(),
        vessel->prolongation().last().degrees_of_freedom());
  }
  new_vessels_.clear();
  LOG(INFO) << "Synchronized the new histories";
}

void Plugin::SynchronizeBubbleHistories() {
  DegreesOfFreedom<Barycentric> const& centre_of_mass =
      current_physics_bubble_->centre_of_mass_trajectory->
          last().degrees_of_freedom();
  for (auto const& pair : current_physics_bubble_->vessels) {
    Vessel* const vessel = pair.first;
    Displacement<Barycentric> const& displacement =
        current_physics_bubble_->displacements_from_centre_of_mass->at(vessel);
    Velocity<Barycentric> const& velocity =
        current_physics_bubble_->velocities_from_centre_of_mass->at(vessel);
    if (vessel->is_synchronized()) {
      vessel->mutable_history()->Append(
          HistoryTime(),
          {centre_of_mass.position + displacement,
           centre_of_mass.velocity + velocity});
    } else {
      vessel->CreateHistoryAndForkProlongation(
          HistoryTime(),
          {centre_of_mass.position + displacement,
           centre_of_mass.velocity + velocity});
      CHECK(new_vessels_.erase(vessel));
    }
  }
}

void Plugin::ResetProlongations() {
  for (auto const& pair : vessels_) {
    std::unique_ptr<Vessel> const& vessel = pair.second;
    vessel->ResetProlongation(HistoryTime());
  }
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    celestial->ResetProlongation(HistoryTime());
  }
  VLOG(1) << "Prolongations have been reset";
}

void Plugin::EvolveProlongationsAndBubble(Instant const& t) {
  NBodySystem<Barycentric>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() + celestials_.size() -
                       NumberOfVesselsInPhysicsBubble() +
                       HavePhysicsBubble() ? 1 : 0);
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (auto const& pair : vessels_) {
    Vessel* const vessel = pair.second.get();
    if (!IsInPhysicsBubble(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  if (HavePhysicsBubble()) {
    trajectories.push_back(
        current_physics_bubble_->centre_of_mass_trajectory.get());
  }
  VLOG(1) << "Evolving prolongations and new histories" << '\n'
          << "from : " << trajectories.front()->last().time() << '\n'
          << "to   : " << t;
  n_body_system_->Integrate(prolongation_integrator_,  // integrator
                            t,                         // tmax
                            Δt_,                       // Δt
                            0,                         // sampling_period
                            true,                      // tmax_is_exact
                            trajectories);             // trajectories
  if (HavePhysicsBubble()) {
    DegreesOfFreedom<Barycentric> const& centre_of_mass =
        current_physics_bubble_->centre_of_mass_trajectory->
            last().degrees_of_freedom();
    for (auto const& pair : current_physics_bubble_->vessels) {
      Vessel* const vessel = pair.first;
      Displacement<Barycentric> const& displacement =
          current_physics_bubble_->
              displacements_from_centre_of_mass->at(vessel);
      Velocity<Barycentric> const& velocity =
          current_physics_bubble_->velocities_from_centre_of_mass->at(vessel);
      vessel->mutable_prolongation()->Append(
          t,
          {centre_of_mass.position + displacement,
           centre_of_mass.velocity + velocity});
    }
  }
}

bool Plugin::IsInPhysicsBubble(Vessel* const vessel) const {
  return HavePhysicsBubble() &&
         current_physics_bubble_->vessels.find(vessel) !=
             current_physics_bubble_->vessels.end();
}

std::size_t Plugin::NumberOfVesselsInPhysicsBubble() const {
  if (HavePhysicsBubble()) {
    return current_physics_bubble_->vessels.size();
  } else {
    return 0;
  }
}

bool Plugin::HavePhysicsBubble() const {
  return current_physics_bubble_ != nullptr;
}

Instant const& Plugin::HistoryTime() const {
  return sun_->history().last().time();
}

Plugin::Plugin(Instant const& initial_time,
               Index const sun_index,
               GravitationalParameter const& sun_gravitational_parameter,
               Angle const& planetarium_rotation)
    : n_body_system_(new NBodySystem<Barycentric>),
      planetarium_rotation_(planetarium_rotation),
      current_time_(initial_time) {
  auto inserted = celestials_.emplace(
      sun_index,
      std::make_unique<Celestial>(
          std::make_unique<MassiveBody>(sun_gravitational_parameter)));
  sun_ = inserted.first->second.get();
  sun_->CreateHistoryAndForkProlongation(
      current_time_,
      {Position<Barycentric>(), Velocity<Barycentric>()});
  history_integrator_.Initialize(history_integrator_.Order5Optimal());
  // NOTE(egg): perhaps a lower order would be appropriate.
  prolongation_integrator_.Initialize(history_integrator_.Order5Optimal());
}

void Plugin::InsertCelestial(
    Index const celestial_index,
    GravitationalParameter const& gravitational_parameter,
    Index const parent_index,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun> const& from_parent_velocity) {
  CHECK(initializing) << "Celestial bodies should be inserted before the end "
                      << "of initialization";
  auto const it = celestials_.find(parent_index);
  CHECK(it != celestials_.end()) << "No body at index " << parent_index;
  Celestial const& parent= *it->second;
  auto const inserted = celestials_.emplace(
      celestial_index,
      std::make_unique<Celestial>(
          std::make_unique<MassiveBody>(gravitational_parameter)));
  CHECK(inserted.second) << "Body already exists at index " << celestial_index;
  LOG(INFO) << "Initial |orbit.pos| for celestial at index " << celestial_index
            << ": " << from_parent_position;
  Displacement<Barycentric> const displacement =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_position));
  LOG(INFO) << "In barycentric coordinates: " << displacement;
  LOG(INFO) << "Initial |orbit.vel| for vessel at index " << celestial_index
            << ": " << from_parent_velocity;
  Velocity<Barycentric> const relative_velocity =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_velocity));
  LOG(INFO) << "In barycentric coordinates: " << relative_velocity;
  Celestial* const celestial = inserted.first->second.get();
  celestial->set_parent(&parent);
  auto const last = parent.history().last();
  celestial->CreateHistoryAndForkProlongation(
      current_time_,
      {last.degrees_of_freedom().position + displacement,
       last.degrees_of_freedom().velocity + relative_velocity});
}

void Plugin::EndInitialization() {
  initializing.Flop();
}

void Plugin::UpdateCelestialHierarchy(Index const celestial_index,
                                      Index const parent_index) const {
  CHECK(!initializing);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  auto const it_parent = celestials_.find(parent_index);
  CHECK(it_parent != celestials_.end()) << "No body at index " << parent_index;
  it->second->set_parent(it_parent->second.get());
}

bool Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                Index const parent_index) {
  CHECK(!initializing);
  auto const it = celestials_.find(parent_index);
  CHECK(it != celestials_.end()) << "No body at index " << parent_index;
  Celestial const& parent = *it->second;
  auto inserted = vessels_.emplace(vessel_guid,
                                   std::make_unique<Vessel>(&parent));
  Vessel* const vessel = inserted.first->second.get();
  kept_.emplace(vessel);
  vessel->set_parent(&parent);
  LOG_IF(INFO, inserted.second) << "Inserted Vessel with GUID " << vessel_guid;
  VLOG(1) << "Parent of vessel with GUID " << vessel_guid <<" is at index "
          << parent_index;
  return inserted.second;
}

void Plugin::SetVesselStateOffset(
    GUID const& vessel_guid,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun> const& from_parent_velocity) {
  CHECK(!initializing);
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel* const vessel = it->second.get();
  CHECK(!vessel->is_initialized())
      << "Vessel with GUID " << vessel_guid << " already has a trajectory";
  LOG(INFO) << "Initial |orbit.pos| for vessel with GUID " << vessel_guid
            << ": " << from_parent_position;
  Displacement<Barycentric> const displacement =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_position));
  LOG(INFO) << "In barycentric coordinates: " << displacement;
  LOG(INFO) << "Initial |orbit.vel| for vessel with GUID " << vessel_guid
            << ": " << from_parent_velocity;
  Velocity<Barycentric> const relative_velocity =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_velocity));
  LOG(INFO) << "In barycentric coordinates: " << relative_velocity;
  auto const last = vessel->parent().history().last();
  vessel->CreateProlongation(
      current_time_,
      {last.degrees_of_freedom().position + displacement,
       last.degrees_of_freedom().velocity + relative_velocity});
  auto const inserted = new_vessels_.emplace(vessel);
  CHECK(inserted.second);
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  CHECK(!initializing);
  CleanUpVessels();
  PreparePhysicsBubble(t);
  if (HistoryTime() + Δt_ < t) {
    // The histories are far enough behind that we can advance them at least one
    // step and reset the prolongations.
    EvolveSynchronizedHistories(t);
    if (!new_vessels_.empty() || HavePhysicsBubble()) {
      SynchronizeNewHistoriesAndBubble();
    }
    ResetProlongations();
  }
  EvolveProlongationsAndBubble(t);
  VLOG(1) << "Time has been advanced" << '\n'
          << "from : " << current_time_ << '\n'
          << "to   : " << t;
  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
}

Displacement<AliceSun> Plugin::VesselDisplacementFromParent(
    GUID const& vessel_guid) const {
  CHECK(!initializing);
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel const& vessel = *it->second;
  CHECK(vessel.is_initialized()) << "Vessel with GUID " << vessel_guid
                                 << " was not given an initial state";
  Displacement<Barycentric> const barycentric_result =
      vessel.prolongation().last().degrees_of_freedom().position -
      vessel.parent().prolongation().last().degrees_of_freedom().position;
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Vessel with GUID " << vessel_guid << " is at parent position + "
          << barycentric_result << " Barycentre (" << result << " AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::VesselParentRelativeVelocity(
    GUID const& vessel_guid) const {
  CHECK(!initializing);
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel const& vessel = *it->second;
  CHECK(vessel.is_initialized()) << "Vessel with GUID " << vessel_guid
                                 << " was not given an initial state";
  Velocity<Barycentric> const barycentric_result =
      vessel.prolongation().last().degrees_of_freedom().velocity -
      vessel.parent().prolongation().last().degrees_of_freedom().velocity;
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Vessel with GUID " << vessel_guid
          << " moves at parent velocity + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

Displacement<AliceSun> Plugin::CelestialDisplacementFromParent(
    Index const celestial_index) const {
  CHECK(!initializing);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  Celestial const& celestial = *it->second;
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  Displacement<Barycentric> const barycentric_result =
      celestial.prolongation().last().degrees_of_freedom().position -
      celestial.parent().prolongation().last().degrees_of_freedom().position;
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " is at parent position + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::CelestialParentRelativeVelocity(
    Index const celestial_index) const {
  CHECK(!initializing);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  Celestial const& celestial = *it->second;
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  Velocity<Barycentric> const barycentric_result =
      celestial.prolongation().last().degrees_of_freedom().velocity -
      celestial.parent().prolongation().last().degrees_of_freedom().velocity;
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " moves at parent velocity + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

RenderedTrajectory<World> Plugin::RenderedVesselTrajectory(
    GUID const& vessel_guid,
    RenderingFrame const& frame,
    Position<World> const& sun_world_position) const {
  CHECK(!initializing);
  auto const to_world =
      AffineMap<Barycentric, World, Length, Rotation>(
          sun_->prolongation().last().degrees_of_freedom().position,
          sun_world_position,
          Rotation<WorldSun, World>::Identity() * PlanetariumRotation());
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end());
  Vessel const& vessel = *(it->second);
  CHECK(vessel.is_initialized());
  VLOG(1) << "Rendering a trajectory for the vessel with GUID " << vessel_guid;
  RenderedTrajectory<World> result;
  if (!vessel.is_synchronized()) {
    // TODO(egg): We render neither unsynchronized histories nor prolongations
    // at the moment.
    VLOG(1) << "Returning an empty trajectory";
    return result;
  }
  DegreesOfFreedom<Barycentric> const* initial_state = nullptr;
  DegreesOfFreedom<Barycentric> const* final_state = nullptr;
  std::unique_ptr<Trajectory<Barycentric>> const apparent_trajectory =
      frame.ApparentTrajectory(vessel.history());
  for (Trajectory<Barycentric>::NativeIterator
           it = apparent_trajectory->first();
       !it.at_end();
       ++it) {
    final_state = &it.degrees_of_freedom();
    if (initial_state != nullptr) {
      result.emplace_back(to_world(initial_state->position),
                          to_world(final_state->position));
    }
    std::swap(final_state, initial_state);
  }
  VLOG(1) << "Returning a " << result.size() << "-segment trajectory";
  return result;
}

std::unique_ptr<BodyCentredNonRotatingFrame>
Plugin::NewBodyCentredNonRotatingFrame(Index const reference_body_index) const {
  auto const it = celestials_.find(reference_body_index);
  CHECK(it != celestials_.end());
  Celestial const& reference_body = *it->second;
  return std::make_unique<BodyCentredNonRotatingFrame>(reference_body);
}

std::unique_ptr<BarycentricRotatingFrame> Plugin::NewBarycentricRotatingFrame(
    Index const primary_index,
    Index const secondary_index) const {
  auto const primary_it = celestials_.find(primary_index);
  CHECK(primary_it != celestials_.end());
  Celestial const& primary = *primary_it->second;
  auto const secondary_it = celestials_.find(secondary_index);
  CHECK(secondary_it != celestials_.end());
  Celestial const& secondary = *secondary_it->second;
  return std::make_unique<BarycentricRotatingFrame>(primary, secondary);
}

Position<World> Plugin::VesselWorldPosition(
    GUID const& vessel_guid,
    Position<World> const& parent_world_position) const {
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end());
  Vessel const& vessel = *(it->second);
  auto const to_world =
      AffineMap<Barycentric, World, Length, Rotation>(
          vessel.parent().prolongation().last().degrees_of_freedom().position,
          parent_world_position,
          Rotation<WorldSun, World>::Identity() * PlanetariumRotation());
  CHECK(vessel.is_initialized()) << "Vessel with GUID " << vessel_guid
                                 << " was not given an initial state";
  return to_world(
      vessel.prolongation().last().degrees_of_freedom().position);
}

Velocity<World> Plugin::VesselWorldVelocity(
      GUID const& vessel_guid,
      Velocity<World> const& parent_world_velocity,
      Time const& parent_rotation_period) const {
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end());
  Vessel const& vessel = *(it->second);
  CHECK(vessel.is_initialized()) << "Vessel with GUID " << vessel_guid
                              << " was not given an initial state";
  Rotation<Barycentric, World> to_world =
      Rotation<WorldSun, World>::Identity() * PlanetariumRotation();
  Velocity<Barycentric> const velocity_relative_to_parent =
      vessel.prolongation().last().degrees_of_freedom().velocity -
      vessel.parent().prolongation().last().degrees_of_freedom().velocity;
  Displacement<Barycentric> const offset_from_parent =
      vessel.prolongation().last().degrees_of_freedom().position -
      vessel.parent().prolongation().last().degrees_of_freedom().position;
  AngularVelocity<Barycentric> const world_frame_angular_velocity =
      AngularVelocity<Barycentric>({0 * Radian / Second,
                                    2 * π * Radian / parent_rotation_period,
                                    0 * Radian / Second});
  return to_world(
      (world_frame_angular_velocity * offset_from_parent) / Radian
          + velocity_relative_to_parent) + parent_world_velocity;
}

void Plugin::AddVesselToNextPhysicsBubble(
    GUID const& vessel_guid,
    std::vector<std::pair<PartID, std::unique_ptr<Part<World>>>> parts) {
  if (next_physics_bubble_ == nullptr) {
    next_physics_bubble_ = std::make_unique<PhysicsBubble>();
  }
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end());
  Vessel* const vessel = it->second.get();
  auto const inserted_vessel =
      next_physics_bubble_->vessels.emplace(vessel,
                                            std::vector<Part<World>* const>());
  CHECK(inserted_vessel.second);
  std::vector<Part<World>* const>* const vessel_parts =
      &inserted_vessel.first->second;
  for (std::pair<PartID, std::unique_ptr<Part<World>>>& id_part : parts) {
    CHECK(inserted_vessel.second);
    auto const inserted_part =
    next_physics_bubble_->parts.insert(std::move(id_part));
    CHECK(inserted_part.second);
    vessel_parts->push_back(inserted_part.first->second.get());
  }
}

Displacement<World> Plugin::BubbleDisplacementOffset(
    Position<World> const& sun_world_position) const {
  CHECK(HavePhysicsBubble());
  return Identity<WorldSun, World>()(PlanetariumRotation()(
             current_physics_bubble_->centre_of_mass_trajectory->
                 last().degrees_of_freedom().position -
             sun_->prolongation().last().degrees_of_freedom().position)) +
         sun_world_position - current_physics_bubble_->centre_of_mass->position;
}

Velocity<World> Plugin::BubbleVelocityOffset(
    Index const reference_body_index) const {
  CHECK(HavePhysicsBubble());
  auto const found = celestials_.find(reference_body_index);
  CHECK(found != celestials_.end());
  Celestial const& reference_body = *found->second;
  return Identity<WorldSun, World>()(PlanetariumRotation()(
             current_physics_bubble_->centre_of_mass_trajectory->
                 last().degrees_of_freedom().velocity -
             reference_body.prolongation().
                 last().degrees_of_freedom().velocity));
}

void Plugin::RestartNextPhysicsBubble() {
  CHECK(next_physics_bubble_ != nullptr);
  std::vector<DegreesOfFreedom<Barycentric>> vessel_degrees_of_freedom;
  vessel_degrees_of_freedom.reserve(next_physics_bubble_->vessels.size());
  std::vector<Mass> vessel_masses;
  vessel_masses.reserve(next_physics_bubble_->vessels.size());
  for (auto const& vessel_parts : next_physics_bubble_->vessels) {
    vessel_degrees_of_freedom.push_back(
        vessel_parts.first->
            prolongation().last().degrees_of_freedom());
    vessel_masses.push_back(Mass());
    for (Part<World> const* const part : vessel_parts.second) {
      vessel_masses.back() += part->mass;
    }
  }
  next_physics_bubble_->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(bubble_body_);
  next_physics_bubble_->centre_of_mass_trajectory->Append(
      current_time_,
      physics::Barycentre(vessel_degrees_of_freedom, vessel_masses));
}

Vector<Acceleration, World> Plugin::IntrinsicAcceleration(
    Instant const& next_time,
    std::vector<std::pair<Part<World>*, Part<World>*>>* const common_parts) {
  CHECK_NOTNULL(common_parts);
  CHECK(common_parts->empty());
  // Most of the time no parts explode.  We reserve accordingly.
  common_parts->reserve(current_physics_bubble_->parts.size());
  Vector<Force, World> weighted_sum;
  Mass total_mass;
  auto it_in_current_parts = current_physics_bubble_->parts.cbegin();
  auto it_in_next_parts = next_physics_bubble_->parts.cbegin();
  while (it_in_current_parts != current_physics_bubble_->parts.end() &&
         it_in_next_parts != next_physics_bubble_->parts.end()) {
    if (it_in_current_parts->first < it_in_next_parts->first) {
      ++it_in_current_parts;
    } else if (it_in_next_parts->first < it_in_current_parts->first) {
      ++it_in_next_parts;
    } else {
      Part<World>* const current_part = it_in_current_parts->second.get();
      Part<World>* const next_part = it_in_next_parts->second.get();
      common_parts->emplace_back(current_part, next_part);
      // TODO(egg): not sure what we actually want to do here.
      Mass const mass = (next_part->mass + current_part->mass) / 2.0;
      weighted_sum += ((next_part->degrees_of_freedom.velocity -
                        current_part->degrees_of_freedom.velocity) /
                           (next_time - current_time_) -
                       current_part->expected_ksp_gravity) * mass;
      total_mass += mass;
      ++it_in_current_parts;
      ++it_in_next_parts;
    }
  }
  return weighted_sum / total_mass;
}

void Plugin::ShiftBubble(
    std::vector<std::pair<Part<World>*,
                          Part<World>*>> const* const common_parts) {
  CHECK_NOTNULL(common_parts);
  CHECK(next_physics_bubble_ != nullptr);
  std::vector<DegreesOfFreedom<World>> current_common_degrees_of_freedom;
  current_common_degrees_of_freedom.reserve(common_parts->size());
  std::vector<Mass> current_common_masses;
  current_common_masses.reserve(common_parts->size());
  std::vector<DegreesOfFreedom<World>> next_common_degrees_of_freedom;
  next_common_degrees_of_freedom.reserve(common_parts->size());
  std::vector<Mass> next_common_masses;
  next_common_masses.reserve(common_parts->size());
  for (auto const& current_next : *common_parts) {
    current_common_degrees_of_freedom.emplace_back(
        current_next.first->degrees_of_freedom);
    current_common_masses.emplace_back(current_next.first->mass);
    next_common_degrees_of_freedom.emplace_back(
        current_next.second->degrees_of_freedom);
    next_common_masses.emplace_back(current_next.second->mass);
  }
  DegreesOfFreedom<World> const current_common_centre_of_mass =
      physics::Barycentre(current_common_degrees_of_freedom,
                          current_common_masses);
  DegreesOfFreedom<World> const next_common_centre_of_mass =
      physics::Barycentre(next_common_degrees_of_freedom, next_common_masses);
  // The change in the position of the overall centre of mass resulting from
  // fixing the centre of mass of the intersection.
  Displacement<World> const position_change =
      (next_physics_bubble_->centre_of_mass->position -
           next_common_centre_of_mass.position) -
      (current_physics_bubble_->centre_of_mass->position -
           current_common_centre_of_mass.position);
  // The change in the velocity of the overall centre of mass resulting from
  // fixing the velocity of the centre of mass of the intersection.
  Velocity<World> const velocity_change =
      (next_physics_bubble_->centre_of_mass->velocity -
           next_common_centre_of_mass.velocity) -
      (current_physics_bubble_->centre_of_mass->velocity -
           current_common_centre_of_mass.velocity);
  DegreesOfFreedom<Barycentric> const& current_centre_of_mass =
      current_physics_bubble_->
          centre_of_mass_trajectory->last().degrees_of_freedom();
  next_physics_bubble_->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(bubble_body_);
  // Using the identity as the map |World| -> |WorldSun| is valid for
  // velocities too since we assume |World| is currently nonrotating, i.e.,
  // it is stationary with respect |to WorldSun|.
  next_physics_bubble_->centre_of_mass_trajectory->Append(
      current_time_,
      {current_centre_of_mass.position +
           PlanetariumRotation().Inverse()(
               Identity<World, WorldSun>()(position_change)),
       current_centre_of_mass.velocity +
           PlanetariumRotation().Inverse()(
               Identity<World, WorldSun>()(velocity_change))});
}

void Plugin::ComputeNextPhysicsBubbleCentreOfMassWorldDegreesOfFreedom() {
  CHECK(next_physics_bubble_ != nullptr);
  std::vector<DegreesOfFreedom<World>> part_degrees_of_freedom;
  part_degrees_of_freedom.reserve(next_physics_bubble_->parts.size());
  std::vector<Mass> part_masses;
  part_masses.reserve(next_physics_bubble_->parts.size());
  for (auto const& id_part : next_physics_bubble_->parts) {
    part_degrees_of_freedom.push_back(id_part.second->degrees_of_freedom);
    part_masses.push_back(id_part.second->mass);
  }
  next_physics_bubble_->centre_of_mass =
      std::make_unique<DegreesOfFreedom<World>>(
          physics::Barycentre(part_degrees_of_freedom, part_masses));
}

void Plugin::ComputeNextPhysicsBubbleVesselOffsets() {
  VLOG(1) << "ComputeNextPhysicsBubbleVesselOffsets";
  CHECK(next_physics_bubble_ != nullptr);
  next_physics_bubble_->displacements_from_centre_of_mass =
      std::make_unique<std::map<Vessel const* const,
                                Displacement<Barycentric>>>();
  next_physics_bubble_->velocities_from_centre_of_mass =
      std::make_unique<std::map<Vessel const* const,
                                Velocity<Barycentric>>>();
  VLOG(1) << "Vessels in next bubble: "
          << next_physics_bubble_->vessels.size();
  for (auto const& vessel_parts : next_physics_bubble_->vessels) {
    VLOG(1) << "Parts in vessel: " << vessel_parts.second.size();
    std::vector<DegreesOfFreedom<World>> part_degrees_of_freedom;
    std::vector<Mass> part_masses;
    part_degrees_of_freedom.reserve(vessel_parts.second.size());
    part_masses.reserve(vessel_parts.second.size());
    for (auto const part : vessel_parts.second) {
      part_degrees_of_freedom.emplace_back(part->degrees_of_freedom);
      part_masses.emplace_back(part->mass);
    }
    DegreesOfFreedom<World> const vessel_degrees_of_freedom =
        physics::Barycentre(part_degrees_of_freedom, part_masses);
    Displacement<Barycentric> const displacement =
        PlanetariumRotation().Inverse()(
            Identity<World, WorldSun>()(
                vessel_degrees_of_freedom.position -
                next_physics_bubble_->centre_of_mass->position));
    Velocity<Barycentric> const velocity =
        PlanetariumRotation().Inverse()(
            Identity<World, WorldSun>()(
                vessel_degrees_of_freedom.velocity -
                next_physics_bubble_->centre_of_mass->velocity));
    next_physics_bubble_->displacements_from_centre_of_mass->emplace(
        vessel_parts.first,
        displacement);
    next_physics_bubble_->velocities_from_centre_of_mass->emplace(
        vessel_parts.first,
        velocity);
  }
}

void Plugin::PreparePhysicsBubble(Instant const& next_time) {
  VLOG(1) << "PreparePhysicsBubble";
  if (next_physics_bubble_ != nullptr) {
    ComputeNextPhysicsBubbleCentreOfMassWorldDegreesOfFreedom();
    ComputeNextPhysicsBubbleVesselOffsets();
    if (current_physics_bubble_ == nullptr) {
      // There was no physics bubble.
      RestartNextPhysicsBubble();
    } else {
      // The IDs of the parts that are both in the current and in the next
      // physics bubble.
      auto const common_parts =
          std::make_unique<
              std::vector<std::pair<Part<World>*, Part<World>*>>>();
      Vector<Acceleration, World> intrinsic_acceleration =
          IntrinsicAcceleration(next_time, common_parts.get());
      if (common_parts->empty()) {
        // The current and next set of parts are disjoint, i.e., the next
        // physics bubble is unrelated to the current one.
        RestartNextPhysicsBubble();
      } else {
        if (common_parts->size() == next_physics_bubble_->parts.size()) {
          // The set of parts has not changed.
          next_physics_bubble_->centre_of_mass_trajectory =
              std::move(current_physics_bubble_->centre_of_mass_trajectory);
          // TODO(egg): we end up dragging some history along here, we probably
          // should not.
        } else {
          // Parts appeared or were removed from the physics bubble, but the
          // intersection is nonempty.  We fix the degrees of freedom of the
          // centre of mass of the intersection, and we use its measured
          // acceleration as the intrinsic acceleration of the |bubble_body_|.
          ShiftBubble(common_parts.get());
        }
        // Correct since |World| is currently nonrotating.
        Vector<Acceleration, Barycentric> barycentric_intrinsic_acceleration =
            PlanetariumRotation().Inverse()(
                Identity<World, WorldSun>()(intrinsic_acceleration));
        LOG(INFO) << "Intrinsic accerelation: "
                  << barycentric_intrinsic_acceleration;
        if (next_physics_bubble_->centre_of_mass_trajectory->
                has_intrinsic_acceleration()) {
          next_physics_bubble_->centre_of_mass_trajectory->
              clear_intrinsic_acceleration();
        }
        // TODO(egg): this makes the intrinsic acceleration a step function.  Might
        // something smoother be better?  We need to be careful not to be one step
        // or half a step in the past though.
        next_physics_bubble_->centre_of_mass_trajectory->
              set_intrinsic_acceleration(
                  [barycentric_intrinsic_acceleration](Instant const& t) {
                    return barycentric_intrinsic_acceleration;
                  });
      }
    }
  }
  current_physics_bubble_ = std::move(next_physics_bubble_);
}

}  // namespace ksp_plugin
}  // namespace principia
