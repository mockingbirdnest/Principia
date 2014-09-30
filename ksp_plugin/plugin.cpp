#include "ksp_plugin/plugin.hpp"

#include <cmath>
#include <string>

#include "geometry/permutation.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::Bivector;
using geometry::Permutation;

Permutation<World, AliceWorld> const kWorldLookingGlass(
    Permutation<World, AliceWorld>::CoordinatePermutation::XZY);

Permutation<WorldSun, AliceSun> const kSunLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

void Plugin::CheckVesselInvariants(
    Vessel const& vessel,
    GUIDToUnownedVessel::iterator const it_in_new_vessels) const {
  CHECK(vessel.history != nullptr) << "Vessel with GUID "
                                   << it_in_new_vessels->first
                                   << " was not given an initial state";
  // TODO(egg): At the moment, if a vessel is inserted when
  // |current_time_ == HistoryTime()| (that only happens before the first call
  // to |AdvanceTime|) its first step is unsynchronised. This is convenient to
  // test code paths, but it means the invariant is GE, rather than GT.
  if (it_in_new_vessels != new_vessels_.end()) {
    CHECK(vessel.prolongation == nullptr);
    CHECK_GE(vessel.history->last_time(), HistoryTime());
  } else {
    CHECK_NOTNULL(vessel.prolongation);
    CHECK_EQ(vessel.history->last_time(), HistoryTime());
  }
}

void Plugin::CleanUpVessels() {
  VLOG(1) << "Vessel cleanup";
  // Remove the vessels which were not updated since last time.
  for (auto& it = vessels_.cbegin(); it != vessels_.cend();) {
    auto const& it_in_new_vessels = new_vessels_.find(it->first);
    Vessel const& vessel = *it->second;
    // While we're going over the vessels, check invariants.
    CheckVesselInvariants(vessel, it_in_new_vessels);
    // Now do the cleanup.
    if (vessel.keep) {
      it->second->keep = false;
      ++it;
    } else {
      LOG(INFO) << "Removing vessel with GUID " << it->first;
      // Since we are going to delete the vessel, we must remove it from
      // |new_vessels| if it's there.
      if (it_in_new_vessels != new_vessels_.end()) {
        LOG(INFO) << "Vessel had not been synchronised";
        new_vessels_.erase(it_in_new_vessels);
      }
      // |std::map::erase| invalidates its parameter so we post-increment.
      vessels_.erase(it++);
    }
  }
}

void Plugin::EvolveSynchronisedHistories(Instant const& t) {
  VLOG(1) << "Starting the evolution of the old histories" << '\n'
          << "from : " << HistoryTime();
  // Integration with a constant step.
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() - new_vessels_.size() +
                       celestials_.size());
  for (auto const& pair : celestials_) {
    trajectories.push_back(pair.second->history.get());
  }
  for (auto const& pair : vessels_) {
    if (pair.second->prolongation != nullptr) {
      trajectories.push_back(pair.second->history.get());
    }
  }
  solar_system_.Integrate(history_integrator_,  // integrator
                          t,                    // tmax
                          kΔt,                  // Δt
                          0,                    // sampling_period
                          false,                // tmax_is_exact
                          trajectories);        // trajectories
  CHECK_GT(HistoryTime(), current_time_);
  VLOG(1) << "Evolved the old histories" << '\n'
          << "to   : " << HistoryTime();
}

void Plugin::SynchroniseNewHistories() {
  VLOG(1) << "Starting the synchronisation of the new histories";
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(celestials_.size() + new_vessels_.size());
  for (auto const& pair : celestials_) {
    trajectories.push_back(pair.second->prolongation);
  }
  for (auto const& pair : new_vessels_) {
    trajectories.push_back(pair.second->history.get());
  }
  solar_system_.Integrate(prolongation_integrator_,  // integrator
                          HistoryTime(),             // tmax
                          kΔt,                       // Δt
                          0,                         // sampling_period
                          true,                      // tmax_is_exact
                          trajectories);             // trajectories
  new_vessels_.clear();
  LOG(INFO) << "Synchronised the new histories";
}

void Plugin::ResetProlongations() {
  for (auto const& pair : vessels_) {
    if (pair.second->prolongation != nullptr) {
      pair.second->history->DeleteFork(&pair.second->prolongation);
    }
    pair.second->prolongation = pair.second->history->Fork(HistoryTime());
  }
  for (auto const& pair : celestials_) {
    pair.second->history->DeleteFork(&pair.second->prolongation);
    pair.second->prolongation = pair.second->history->Fork(HistoryTime());
  }
  VLOG(1) << "Prolongations have been reset";
}

void Plugin::EvolveProlongationsAndUnsynchronisedHistories(Instant const& t) {
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() + celestials_.size());
  for (auto const& pair : celestials_) {
    trajectories.push_back(pair.second->prolongation);
  }
  for (auto const& pair : new_vessels_) {
    trajectories.push_back(pair.second->history.get());
  }
  for (auto const& pair : vessels_) {
    if (pair.second->prolongation != nullptr) {
      trajectories.push_back(pair.second->prolongation);
    }
  }
  VLOG(1) << "Evolving prolongations and new histories" << '\n'
          << "from : " << trajectories.front()->last_time() << '\n'
          << "to   : " << t;
  solar_system_.Integrate(prolongation_integrator_,  // integrator
                          t,                         // tmax
                          kΔt,                       // Δt
                          0,                         // sampling_period
                          true,                      // tmax_is_exact
                          trajectories);             // trajectories
}

Instant const& Plugin::HistoryTime() const {
  return sun_->history->last_time();
}

Rotation<Barycentre, WorldSun> Plugin::PlanetariumRotation() const {
  return Rotation<Barycentre, WorldSun>(
      planetarium_rotation_,
      Bivector<double, Barycentre>({0, 1, 0}));
}

Plugin::Plugin(Instant const& initial_time,
               Index const sun_index,
               GravitationalParameter const& sun_gravitational_parameter,
               Angle const& planetarium_rotation)
    : current_time_(initial_time),
      planetarium_rotation_(planetarium_rotation) {
  auto inserted = celestials_.insert(
      {sun_index, std::make_unique<Celestial>(sun_gravitational_parameter)});
  sun_ = inserted.first->second.get();
  sun_->history = std::make_unique<Trajectory<Barycentre>>(*sun_->body);
  sun_->history->Append(current_time_,
                        {Position<Barycentre>(), Velocity<Barycentre>()});
  sun_->prolongation = sun_->history->Fork(current_time_);
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
  CHECK(initialising) << "Celestial bodies should be inserted before the first "
                      << "call to |AdvanceTime|";
  auto const it = celestials_.find(parent_index);
  CHECK(it != celestials_.end()) << "No body at index " << parent_index;
  Celestial const& parent= *it->second;
  auto const inserted = celestials_.insert(
      {celestial_index, std::make_unique<Celestial>(gravitational_parameter)});
  CHECK(inserted.second) << "Body already exists at index " << celestial_index;
  LOG(INFO) << "Initial |orbit.pos| for celestial at index " << celestial_index
            << ": " << from_parent_position;
  Displacement<Barycentre> const displacement =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_position));
  LOG(INFO) << "In barycentric coordinates: " << displacement;
  LOG(INFO) << "Initial |orbit.vel| for vessel at index " << celestial_index
            << ": " << from_parent_velocity;
  Velocity<Barycentre> const relative_velocity =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_velocity));
  LOG(INFO) << "In barycentric coordinates: " << relative_velocity;
  Celestial* const celestial = inserted.first->second.get();
  celestial->parent = &parent;
  celestial->history =
      std::make_unique<Trajectory<Barycentre>>(*celestial->body);
  celestial->history->Append(
      current_time_,
      {parent.history->last_position() + displacement,
       parent.history->last_velocity() + relative_velocity});
  celestial->prolongation = celestial->history->Fork(current_time_);
}

void Plugin::EndInitialisation() {
  initialising.Flop();
}

void Plugin::UpdateCelestialHierarchy(Index const celestial_index,
                                      Index const parent_index) const {
  CHECK(!initialising);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  auto const it_parent = celestials_.find(parent_index);
  CHECK(it_parent != celestials_.end()) << "No body at index " << parent_index;
  it->second->parent = it_parent->second.get();
}

bool Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                Index const parent_index) {
  CHECK(!initialising);
  auto const it = celestials_.find(parent_index);
  CHECK(it != celestials_.end()) << "No body at index " << parent_index;
  Celestial const& parent = *it->second;
  auto inserted = vessels_.insert({vessel_guid,
                                   std::make_unique<Vessel>(&parent)});
  Vessel* const vessel = inserted.first->second.get();
  vessel->keep = true;
  vessel->parent = &parent;
  LOG_IF(INFO, inserted.second) << "Inserted Vessel with GUID " << vessel_guid;
  VLOG(1) << "Parent of vessel with GUID " << vessel_guid <<" is at index "
          << parent_index;
  return inserted.second;
}

void Plugin::SetVesselStateOffset(
    GUID const& vessel_guid,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun> const& from_parent_velocity) {
  CHECK(!initialising);
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel* const vessel = it->second.get();
  CHECK(vessel->history == nullptr)
      << "Vessel with GUID " << vessel_guid << " already has a trajectory";
  LOG(INFO) << "Initial |orbit.pos| for vessel with GUID " << vessel_guid
            << ": " << from_parent_position;
  Displacement<Barycentre> const displacement =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_position));
  LOG(INFO) << "In barycentric coordinates: " << displacement;
  LOG(INFO) << "Initial |orbit.vel| for vessel with GUID " << vessel_guid
            << ": " << from_parent_velocity;
  Velocity<Barycentre> const relative_velocity =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_velocity));
  LOG(INFO) << "In barycentric coordinates: " << relative_velocity;
  vessel->history = std::make_unique<Trajectory<Barycentre>>(*vessel->body);
  vessel->history->Append(
      current_time_,
      {vessel->parent->history->last_position() + displacement,
       vessel->parent->history->last_velocity() + relative_velocity});
  auto const inserted = new_vessels_.insert({vessel_guid, vessel});
  CHECK(inserted.second);
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  CHECK(!initialising);
  CleanUpVessels();
  if (HistoryTime() + kΔt < t) {
    // The histories are far enough behind that we can advance them at least one
    // step and reset the prolongations.
    EvolveSynchronisedHistories(t);
    if (!new_vessels_.empty()) {
      SynchroniseNewHistories();
    }
    ResetProlongations();
  }
  EvolveProlongationsAndUnsynchronisedHistories(t);
  VLOG(1) << "Time has been advanced" << '\n'
          << "from : " << current_time_ << '\n'
          << "to   : " << t;
  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
}

Displacement<AliceSun> Plugin::VesselDisplacementFromParent(
    GUID const& vessel_guid) const {
  CHECK(!initialising);
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel const& vessel = *it->second;
  CHECK(vessel.history != nullptr) << "Vessel with GUID " << vessel_guid
                                   << " was not given an initial state";
  Displacement<Barycentre> const barycentric_result =
      (vessel.prolongation == nullptr ?
           vessel.history->last_position() :
           vessel.prolongation->last_position()) -
      vessel.parent->prolongation->last_position();
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Vessel with GUID " << vessel_guid << " is at parent position + "
          << barycentric_result << " Barycentre (" << result << " AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::VesselParentRelativeVelocity(
    GUID const& vessel_guid) const {
  CHECK(!initialising);
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel const& vessel = *it->second;
  CHECK(vessel.history != nullptr) << "Vessel with GUID " << vessel_guid
                                   << " was not given an initial state";
  Velocity<Barycentre> const barycentric_result =
      (vessel.prolongation == nullptr ?
           vessel.history->last_velocity() :
           vessel.prolongation->last_velocity()) -
       vessel.parent->prolongation->last_velocity();
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Vessel with GUID " << vessel_guid
          << " moves at parent velocity + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

Displacement<AliceSun> Plugin::CelestialDisplacementFromParent(
    Index const celestial_index) const {
  CHECK(!initialising);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  Celestial const& celestial = *it->second;
  CHECK(celestial.parent != nullptr)
      << "Body at index " << celestial_index << " is the sun";
  Displacement<Barycentre> const barycentric_result =
      celestial.prolongation->last_position() -
      celestial.parent->prolongation->last_position();
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " is at parent position + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::CelestialParentRelativeVelocity(
    Index const celestial_index) const {
  CHECK(!initialising);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  Celestial const& celestial = *it->second;
  CHECK(celestial.parent != nullptr)
      << "Body at index " << celestial_index << " is the sun";
  Velocity<Barycentre> const barycentric_result =
      celestial.prolongation->last_velocity() -
      celestial.parent->prolongation->last_velocity();
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " moves at parent velocity + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

}  // namespace ksp_plugin
}  // namespace principia
