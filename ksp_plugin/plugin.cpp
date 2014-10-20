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

Permutation<World, AliceWorld> const kWorldLookingGlass(
    Permutation<World, AliceWorld>::CoordinatePermutation::XZY);

Permutation<WorldSun, AliceSun> const kSunLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

Rotation<Barycentre, WorldSun> Plugin::PlanetariumRotation() const {
  return Rotation<Barycentre, WorldSun>(
      planetarium_rotation_,
      Bivector<double, Barycentre>({0, 1, 0}));
}

void Plugin::CheckVesselInvariants(
    Vessel<Barycentre> const& vessel,
    GUIDToUnownedVessel::iterator const it_in_new_vessels) const {
  CHECK(vessel.has_history()) << "Vessel with GUID " << it_in_new_vessels->first
                              << " was not given an initial state";
  // TODO(egg): At the moment, if a vessel is inserted when
  // |current_time_ == HistoryTime()| (that only happens before the first call
  // to |AdvanceTime|) its first step is unsynchronized. This is convenient to
  // test code paths, but it means the invariant is GE, rather than GT.
  if (it_in_new_vessels != new_vessels_.end()) {
    CHECK(!vessel.has_prolongation());
    CHECK_GE(vessel.history().last_time(), HistoryTime());
  } else {
    CHECK(vessel.has_prolongation());
    CHECK_EQ(vessel.history().last_time(), HistoryTime());
  }
}

void Plugin::CleanUpVessels() {
  VLOG(1) << "Vessel cleanup";
  // Remove the vessels which were not updated since last time.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    auto const& it_in_new_vessels = new_vessels_.find(it->first);
    Vessel<Barycentre> const* vessel = it->second.get();
    // While we're going over the vessels, check invariants.
    CheckVesselInvariants(*vessel, it_in_new_vessels);
    // Now do the cleanup.
    if (kept_.erase(vessel)) {
      ++it;
    } else {
      LOG(INFO) << "Removing vessel with GUID " << it->first;
      // Since we are going to delete the vessel, we must remove it from
      // |new_vessels| if it's there.
      if (it_in_new_vessels != new_vessels_.end()) {
        LOG(INFO) << "Vessel had not been synchronized";
        new_vessels_.erase(it_in_new_vessels);
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
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() - new_vessels_.size() +
                       celestials_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial<Barycentre>> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_history());
  }
  for (auto const& pair : vessels_) {
    std::unique_ptr<Vessel<Barycentre>> const& vessel = pair.second;
    if (vessel->has_prolongation()) {
      trajectories.push_back(vessel->mutable_history());
    }
  }
  n_body_system_->Integrate(history_integrator_,  // integrator
                            t,                    // tmax
                            Δt_,                  // Δt
                            0,                    // sampling_period
                            false,                // tmax_is_exact
                            trajectories);        // trajectories
  CHECK_GT(HistoryTime(), current_time_);
  VLOG(1) << "Evolved the old histories" << '\n'
          << "to   : " << HistoryTime();
}

void Plugin::SynchronizeNewHistories() {
  VLOG(1) << "Starting the synchronization of the new histories";
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(celestials_.size() + new_vessels_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial<Barycentre>> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (auto const& pair : new_vessels_) {
    Vessel<Barycentre>* vessel = pair.second;
    trajectories.push_back(vessel->mutable_history());
  }
  n_body_system_->Integrate(prolongation_integrator_,  // integrator
                            HistoryTime(),             // tmax
                            Δt_,                       // Δt
                            0,                         // sampling_period
                            true,                      // tmax_is_exact
                            trajectories);             // trajectories
  new_vessels_.clear();
  LOG(INFO) << "Synchronized the new histories";
}

void Plugin::ResetProlongations() {
  for (auto const& pair : vessels_) {
    std::unique_ptr<Vessel<Barycentre>> const& vessel = pair.second;
    vessel->ResetProlongation(HistoryTime());
  }
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial<Barycentre>> const& celestial = pair.second;
    celestial->ResetProlongation(HistoryTime());
  }
  VLOG(1) << "Prolongations have been reset";
}

void Plugin::EvolveProlongationsAndUnsynchronizedHistories(Instant const& t) {
  NBodySystem<Barycentre>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() + celestials_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial<Barycentre>> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (auto const& pair : new_vessels_) {
    Vessel<Barycentre>* vessel = pair.second;
    trajectories.push_back(vessel->mutable_history());
  }
  for (auto const& pair : vessels_) {
    std::unique_ptr<Vessel<Barycentre>> const& vessel = pair.second;
    if (vessel->has_prolongation()) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  VLOG(1) << "Evolving prolongations and new histories" << '\n'
          << "from : " << trajectories.front()->last_time() << '\n'
          << "to   : " << t;
  n_body_system_->Integrate(prolongation_integrator_,  // integrator
                            t,                         // tmax
                            Δt_,                       // Δt
                            0,                         // sampling_period
                            true,                      // tmax_is_exact
                            trajectories);             // trajectories
}

Instant const& Plugin::HistoryTime() const {
  return sun_->history().last_time();
}

Plugin::Plugin(Instant const& initial_time,
               Index const sun_index,
               GravitationalParameter const& sun_gravitational_parameter,
               Angle const& planetarium_rotation)
    : n_body_system_(new NBodySystem<Barycentre>),
      planetarium_rotation_(planetarium_rotation),
      current_time_(initial_time) {
  auto inserted = celestials_.insert(
      {sun_index,
       std::make_unique<Celestial<Barycentre>>(sun_gravitational_parameter)});
  sun_ = inserted.first->second.get();
  sun_->AppendAndForkProlongation(
      current_time_,
      {Position<Barycentre>(), Velocity<Barycentre>()});
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
  Celestial<Barycentre> const& parent= *it->second;
  auto const inserted = celestials_.insert(
      {celestial_index,
       std::make_unique<Celestial<Barycentre>>(gravitational_parameter)});
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
  Celestial<Barycentre>* const celestial = inserted.first->second.get();
  celestial->set_parent(&parent);
  celestial->AppendAndForkProlongation(
      current_time_,
      {parent.history().last_position() + displacement,
       parent.history().last_velocity() + relative_velocity});
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
  Celestial<Barycentre> const& parent = *it->second;
  auto inserted = vessels_.insert(
      {vessel_guid,
       std::make_unique<Vessel<Barycentre>>(&parent)});
  Vessel<Barycentre>* const vessel = inserted.first->second.get();
  kept_.insert(vessel);
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
  Vessel<Barycentre>* const vessel = it->second.get();
  CHECK(!vessel->has_history())
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
  vessel->Append(
      current_time_,
      {vessel->parent().history().last_position() + displacement,
       vessel->parent().history().last_velocity() + relative_velocity});
  auto const inserted = new_vessels_.insert({vessel_guid, vessel});
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
  EvolveProlongationsAndUnsynchronizedHistories(t);
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
  Vessel<Barycentre> const& vessel = *it->second;
  CHECK(vessel.has_history()) << "Vessel with GUID " << vessel_guid
                              << " was not given an initial state";
  Displacement<Barycentre> const barycentric_result =
      vessel.prolongation_or_history().last_position() -
      vessel.parent().prolongation().last_position();
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
  Vessel<Barycentre> const& vessel = *it->second;
  CHECK(vessel.has_history()) << "Vessel with GUID " << vessel_guid
                              << " was not given an initial state";
  Velocity<Barycentre> const barycentric_result =
      vessel.prolongation_or_history().last_velocity() -
      vessel.parent().prolongation().last_velocity();
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
  Celestial<Barycentre> const& celestial = *it->second;
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  Displacement<Barycentre> const barycentric_result =
      celestial.prolongation().last_position() -
      celestial.parent().prolongation().last_position();
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
  Celestial<Barycentre> const& celestial = *it->second;
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  Velocity<Barycentre> const barycentric_result =
      celestial.prolongation().last_velocity() -
      celestial.parent().prolongation().last_velocity();
  Velocity<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " moves at parent velocity + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

// TODO(egg): would things be faster if we computed a polygon ourselves and had
// Vectrosity render it, rather than telling it to render a spline?
RenderedTrajectory<World> Plugin::RenderedVesselTrajectory(
    GUID const& vessel_guid,
    Instant const& lower_bound,
    RenderingFrame const& frame,
    Position<World> const& sun_world_position) const {
  auto const to_world = AffineMap<Barycentre, World, Length, Rotation>(
        sun_->prolongation().last_position(),
        sun_world_position,
        Rotation<WorldSun, World>::Identity() * PlanetariumRotation());
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end());
  Vessel<Barycentre> const& vessel = *(it->second);
  CHECK(vessel.has_history());
  RenderedTrajectory<World> result;
  if (!vessel.has_prolongation()) {
    // TODO(egg): We render neither unsynchronized histories nor prolongations
    // at the moment.
    return result;
  }
  // Initial and final time and state for the Bézier segment being computed.
  Instant const* initial_time = nullptr;
  Instant const* final_time = nullptr;
  DegreesOfFreedom<Barycentre> const* initial_state = nullptr;
  DegreesOfFreedom<Barycentre> const* final_state = nullptr;
  CubicBézierCurve<Barycentre> barycentric_bézier_points;
  std::unique_ptr<Trajectory<Barycentre>> const apparent_trajectory =
      frame.ApparentTrajectory(vessel.history());
  for (auto const& pair : apparent_trajectory->timeline()) {
    final_time = &pair.first;
    final_state = &pair.second;
    if (initial_state != nullptr && final_time != nullptr) {
      Time const δt = *final_time - *initial_time;
      barycentric_bézier_points[0] = initial_state->position;
      barycentric_bézier_points[3] = final_state->position;
      // TODO(egg): * (1.0 / 3.0) would be faster. Would it be measurably so?
      barycentric_bézier_points[1] =
          barycentric_bézier_points[0] + initial_state->velocity * δt / 3.0;
      barycentric_bézier_points[2] =
          barycentric_bézier_points[3] - final_state->velocity * δt / 3.0;
      result.push_back({to_world(barycentric_bézier_points[0]),
                        to_world(barycentric_bézier_points[1]),
                        to_world(barycentric_bézier_points[2]),
                        to_world(barycentric_bézier_points[3])});
    }
    std::swap(final_time, initial_time);
    std::swap(final_state, initial_state);
  }
  return result;
}

}  // namespace ksp_plugin
}  // namespace principia
