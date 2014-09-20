#include "ksp_plugin/plugin.hpp"

#include <string>

#include "geometry/permutation.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::Bivector;
using geometry::Permutation;

Permutation<WorldSun, AliceSun> const kWorldLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

Permutation<WorldSun, AliceSun> const kSunLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

Rotation<Barycentre, WorldSun> Plugin::PlanetariumRotation() {
  return Rotation<Barycentre, WorldSun>(
      planetarium_rotation_,
      Bivector<double, Barycentre>({0, 1, 0}));
}

Plugin::Plugin(Instant const& initial_time, int const sun_index,
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
    int const index,
    GravitationalParameter const& gravitational_parameter,
    int const parent,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun>  const& from_parent_velocity) {
  CHECK(celestials_.find(parent) != celestials_.end()) << "No body at index "
      << parent;
  Celestial* const parent_body = celestials_[parent].get();
  auto const inserted = celestials_.insert(
      {index, std::make_unique<Celestial>(gravitational_parameter)});
  CHECK(inserted.second) << "Body already exists at index " << index;
  LOG(INFO) << "Initial |orbit.pos| for celestial at index " << index << ": "
            << from_parent_position;
  Displacement<Barycentre> const displacement =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_position));
  LOG(INFO) << "In barycentric coordinates: " << displacement;
  LOG(INFO) << "Initial |orbit.vel| for vessel at index " << index << ": "
            << from_parent_velocity;
  Velocity<Barycentre> const relative_velocity =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_velocity));
  LOG(INFO) << "In barycentric coordinates: " << relative_velocity;
  Celestial* const celestial = inserted.first->second.get();
  celestial->parent = parent_body;
  celestial->history =
      std::make_unique<Trajectory<Barycentre>>(*celestial->body);
  celestial->history->Append(
      current_time_,
      {parent_body->history->last_position() + displacement,
       parent_body->history->last_velocity() + relative_velocity});
}

void Plugin::UpdateCelestialHierarchy(int const index, int const parent) {
  CHECK(celestials_.find(index) != celestials_.end()) <<
      "No body at index " << index;
  CHECK(celestials_.find(parent) != celestials_.end()) <<
      "No body at index " << parent;
  celestials_[index]->parent = celestials_[parent].get();
}

bool Plugin::InsertOrKeepVessel(std::string guid, int const parent) {
  CHECK(celestials_.find(parent) != celestials_.end()) << "No body at index "
      << parent;
  auto inserted = vessels_.insert(
      {guid, std::make_unique<Vessel>(celestials_[parent].get())});
  Vessel* const vessel = inserted.first->second.get();
  vessel->keep = true;
  vessel->parent = celestials_[parent].get();
  LOG(INFO) << (inserted.second ? "Inserted Vessel " : "Vessel ")
            << "with GUID " << guid <<": parent index " << parent;
  return inserted.second;
}

void Plugin::SetVesselStateOffset(
    std::string const& guid,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun> const& from_parent_velocity) {
  CHECK(vessels_.find(guid) != vessels_.end()) << "No vessel with GUID "
      << guid;
  Vessel* const vessel = vessels_[guid].get();
  CHECK(vessel->history == nullptr) << "Vessel with GUID " << guid
      << " already has a trajectory";
  LOG(INFO) << "Initial |orbit.pos| for vessel with GUID " << guid << ": "
            << from_parent_position;
  Displacement<Barycentre> const displacement =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_position));
  LOG(INFO) << "In barycentric coordinates: " << displacement;
  LOG(INFO) << "Initial |orbit.vel| for vessel with GUID " << guid << ": "
            << from_parent_velocity;
  Velocity<Barycentre> const relative_velocity =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent_velocity));
  LOG(INFO) << "In barycentric coordinates: " << relative_velocity;
  vessel->history = std::make_unique<Trajectory<Barycentre>>(*vessel->body);
  vessel->history->Append(
      current_time_,
      {vessel->parent->history->last_position() + displacement,
       vessel->parent->history->last_velocity() + relative_velocity});
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    if (!it->second->keep) {
      // |std::map::erase| invalidates its parameter so we post-increment.
      vessels_.erase(it++);
    } else {
      CHECK(it->second->history != nullptr) << "Vessel with GUID " << it->first
          << "was not given an initial state";
      it->second->keep = false;
      ++it;
    }
  }
  // We make sure the step size divides the interval and is at least as small
  // as |kΔt|.
  Time const duration = t - current_time_;
  Time const Δt = duration / std::ceil(duration / kΔt);
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() + celestials_.size());
  for (auto const& pair : celestials_) {
    trajectories.push_back(pair.second->history.get());
  }
  for (auto const& pair : vessels_) {
    trajectories.push_back(pair.second->history.get());
  }
  solar_system_.Integrate(integrator_, t, Δt, 0, trajectories);
  if (trajectories.front()->last_time() != t) {
    // TODO(egg): This is bound to happen given how the integrator works now...
    LOG(ERROR) << "Integration went "
        << trajectories.front()->last_time() - t << " too far.";
  }
  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
}

Displacement<AliceSun> Plugin::VesselDisplacementFromParent(
    std::string const& guid) {
  CHECK(vessels_.find(guid) != vessels_.end()) << "No vessel with GUID "
      << guid;
  CHECK(vessels_[guid]->history != nullptr) << "Vessel with GUID " << guid
      << "was not given an initial state";
  Displacement<Barycentre> const barycentric_result =
      vessels_[guid]->history->last_position() -
      vessels_[guid]->parent->history->last_position();
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  LOG(INFO) << "Vessel with GUID " << guid << " is at parent position + "
            << barycentric_result << " (Barycentre)";
  LOG(INFO) << result << " (AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::VesselParentRelativeVelocity(
    std::string const& guid) {
  CHECK(vessels_.find(guid) != vessels_.end()) << "No vessel with GUID "
      << guid;
  CHECK(vessels_[guid]->history != nullptr) << "Vessel with GUID " << guid
      << "was not given an initial state";
  Velocity<Barycentre> const barycentric_result =
      vessels_[guid]->history->last_velocity() -
      vessels_[guid]->parent->history->last_velocity();
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  LOG(INFO) << "Vessel with GUID " << guid << " moves at parent velocity + "
            << barycentric_result << " (Barycentre)";
  LOG(INFO) << result << " (AliceSun)";
  return result;
}

Displacement<AliceSun> Plugin::CelestialDisplacementFromParent(int const index) {
  CHECK(celestials_.find(index) != celestials_.end()) << "No body at index "
      << index;
  CHECK(celestials_[index]->parent != nullptr) << "Body at index " << index
      << " is the sun";
  Displacement<Barycentre> const barycentric_result =
      celestials_[index]->history->last_position() -
      celestials_[index]->parent->history->last_position();
  Displacement<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  LOG(INFO) << "Celestial at index " << index << " is at parent position + "
            << barycentric_result << " (Barycentre)";
  LOG(INFO) << result << " (AliceSun)";
  return result;
}

Velocity<AliceSun> Plugin::CelestialParentRelativeVelocity(int const index) {
  CHECK(celestials_.find(index) != celestials_.end()) << "No body at index "
      << index;
  CHECK(celestials_[index]->parent != nullptr) << "Body at index " << index
      << " is the sun";
  Velocity<Barycentre> const barycentric_result =
      celestials_[index]->history->last_velocity() -
      celestials_[index]->parent->history->last_velocity();
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  LOG(INFO) << "Celestial at index " << index << " moves at parent velocity + "
            << barycentric_result << " (Barycentre)";
  LOG(INFO) << result << " (AliceSun)";
  return result;
}

}  // namespace ksp_plugin
}  // namespace principia
