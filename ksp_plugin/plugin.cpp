#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "geometry/permutation.hpp"
#include "geometry/affine_map.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::AffineMap;
using geometry::Bivector;
using geometry::Permutation;
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
  trajectories.reserve(vessels_.size() - new_vessels_.size() +
                       celestials_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_history());
  }
  for (auto const& pair : vessels_) {
    std::unique_ptr<Vessel> const& vessel = pair.second;
    if (vessel->is_synchronized()) {
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

void Plugin::SynchronizeNewHistories() {
  VLOG(1) << "Starting the synchronization of the new histories";
  NBodySystem<Barycentric>::Trajectories trajectories;
  trajectories.reserve(celestials_.size() + new_vessels_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (Vessel* const vessel : new_vessels_) {
    trajectories.push_back(vessel->mutable_prolongation());
  }
  n_body_system_->Integrate(prolongation_integrator_,  // integrator
                            HistoryTime(),             // tmax
                            Δt_,                       // Δt
                            0,                         // sampling_period
                            true,                      // tmax_is_exact
                            trajectories);             // trajectories
  for (Vessel* const vessel : new_vessels_) {
    vessel->CreateHistoryAndForkProlongation(
        HistoryTime(),
        vessel->prolongation().last().degrees_of_freedom());
  }
  new_vessels_.clear();
  LOG(INFO) << "Synchronized the new histories";
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

void Plugin::EvolveProlongations(Instant const& t) {
  NBodySystem<Barycentric>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() + celestials_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (auto const& pair : vessels_) {
    std::unique_ptr<Vessel> const& vessel = pair.second;
    trajectories.push_back(vessel->mutable_prolongation());
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
  if (HistoryTime() + Δt_ < t) {
    // The histories are far enough behind that we can advance them at least one
    // step and reset the prolongations.
    EvolveSynchronizedHistories(t);
    if (!new_vessels_.empty()) {
      SynchronizeNewHistories();
    }
    ResetProlongations();
  }
  EvolveProlongations(t);
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

}  // namespace ksp_plugin
}  // namespace principia
