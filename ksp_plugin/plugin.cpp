#include "ksp_plugin/plugin.hpp"

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
  integrator_.Initialize(integrator_.Order5Optimal());
}

void Plugin::InsertCelestial(
    Index const celestial_index,
    GravitationalParameter const& gravitational_parameter,
    Index const parent_index,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun> const& from_parent_velocity) {
  CHECK(initialising) << "Celestial bodies should be inserted before the first"
                      << " call to |AdvanceTime|";
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
  celestial->rendering_extension = celestial->history->Fork(current_time_);
}

void Plugin::UpdateCelestialHierarchy(Index const celestial_index,
                                      Index const parent_index) const {
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  auto const it_parent = celestials_.find(parent_index);
  CHECK(it_parent != celestials_.end()) << "No body at index " << parent_index;
  it->second->parent = it_parent->second.get();
}

bool Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                Index const parent_index) {
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
  new_vessels_.push_back(vessel);
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  initialising = false;
  // Remove the vessels which were not updated since last time.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    if (!it->second->keep) {
      // |std::map::erase| invalidates its parameter so we post-increment.
      vessels_.erase(it++);
    } else {
      CHECK(it->second->history != nullptr)
          << "Vessel with GUID " << it->first
          << " was not given an initial state";
      it->second->keep = false;
      ++it;
    }
  }
  // Whether there the histories are far enough behind that we can advance them
  // at least one step and reset the rendering extensions.
  // NOTE(egg): 2 is way too large, this should be 1 + ε.
  bool reset_rendering_extensions = history_time_ + 2 * kΔt > t;
  // Integration with a small (non-constant) step.
  if (!reset_rendering_extensions ||  // Extend rendering extensions.
      !new_vessels_.empty()) {        // Catch up histories.
    NBodySystem<Barycentre>::Trajectories trajectories;
    for (auto const& pair : celestials_) {
      trajectories.push_back(pair.second->rendering_extension);
    }
    for (auto const& it : new_vessels_) {
      trajectories.push_back(it->history.get());
    }
    if (!reset_rendering_extensions) {
      for (auto const& pair : vessels_) {
        if (pair.second->rendering_extension != nullptr) {
          trajectories.push_back(pair.second->rendering_extension);
        }
      }
    }
  }
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() + celestials_.size());
  for (auto const& pair : celestials_) {
    trajectories.push_back(pair.second->history.get());
    if (reset_rendering_extension) {
      pair.second->history->DeleteFork(&pair.second->rendering_extension);
    }
  }
  for (auto const& pair : vessels_) {
    trajectories.push_back(pair.second->history.get());
    if (reset_rendering_extension) {
      pair.second->history->DeleteFork(&pair.second->rendering_extension);
    }
  }

  // If histories are far enough behind that we can advance them at least one
  // step, we restart the rendering extensions at the end of the prolonged
  // histories.
  bool const reset_rendering_extension = history_time_ + kΔt < t;
  solar_system_.Integrate(integrator_, t, kΔt, 0, trajectories);
  current_time_ = trajectories.front()->last_time();
  planetarium_rotation_ = planetarium_rotation;
}

Displacement<AliceSun> Plugin::VesselDisplacementFromParent(
    GUID const& vessel_guid) const {
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel const& vessel = *it->second;
  CHECK(vessel.history != nullptr) << "Vessel with GUID " << vessel_guid
                                   << " was not given an initial state";
  Displacement<Barycentre> const barycentric_result =
      vessel.history->last_position() - vessel.parent->history->last_position();
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Vessel with GUID " << vessel_guid << " is at parent position + "
          << barycentric_result << " Barycentre (" << result << " AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::VesselParentRelativeVelocity(
    GUID const& vessel_guid) const {
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  Vessel const& vessel = *it->second;
  CHECK(vessel.history != nullptr) << "Vessel with GUID " << vessel_guid
                                   << " was not given an initial state";
  Velocity<Barycentre> const barycentric_result =
      vessel.history->last_velocity() - vessel.parent->history->last_velocity();
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Vessel with GUID " << vessel_guid
          << " moves at parent velocity + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

Displacement<AliceSun> Plugin::CelestialDisplacementFromParent(
    Index const celestial_index) const {
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  Celestial const& celestial = *it->second;
  CHECK(celestial.parent != nullptr)
      << "Body at index " << celestial_index << " is the sun";
  Displacement<Barycentre> const barycentric_result =
      celestial.history->last_position() -
      celestial.parent->history->last_position();
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " is at parent position + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::CelestialParentRelativeVelocity(
    Index const celestial_index) const {
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  Celestial const& celestial = *it->second;
  CHECK(celestial.parent != nullptr)
      << "Body at index " << celestial_index << " is the sun";
  Velocity<Barycentre> const barycentric_result =
      celestial.history->last_velocity() -
      celestial.parent->history->last_velocity();
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " moves at parent velocity + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

}  // namespace ksp_plugin
}  // namespace principia
