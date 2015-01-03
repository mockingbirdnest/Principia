#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "base/unique_ptr_logging.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/identity.hpp"
#include "geometry/permutation.hpp"
#include "geometry/affine_map.hpp"
#include "glog/logging.h"
#include "glog/stl_logging.h"

namespace principia {
namespace ksp_plugin {

using geometry::AffineMap;
using geometry::BarycentreCalculator;
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

std::unique_ptr<Vessel> const& Plugin::find_vessel_by_guid_or_die(
    GUID const& vessel_guid) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel_guid);
  auto const it = vessels_.find(vessel_guid);
  CHECK(it != vessels_.end()) << "No vessel with GUID " << vessel_guid;
  VLOG_AND_RETURN(1, it->second);
}

bool Plugin::has_dirty_vessels() const {
  return !dirty_vessels_.empty();
}

bool Plugin::has_unsynchronized_vessels() const {
  return !unsynchronized_vessels_.empty();
}

bool Plugin::is_dirty(Vessel* const vessel) const {
  return dirty_vessels_.count(vessel) > 0;
}

// The map between the vector spaces of |Barycentric| and |WorldSun| at
// |current_time_|.
Rotation<Barycentric, WorldSun> Plugin::PlanetariumRotation() const {
  return Rotation<Barycentric, WorldSun>(
      planetarium_rotation_,
      Bivector<double, Barycentric>({0, 1, 0}));
}

void Plugin::CheckVesselInvariants(
    GUIDToOwnedVessel::const_iterator const it) const {
  std::unique_ptr<Vessel> const& vessel = it->second;
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

void Plugin::CleanUpVessels() {
  VLOG(1) <<  __FUNCTION__;
  // Remove the vessels which were not updated since last time.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    // While we're going over the vessels, check invariants.
    CheckVesselInvariants(it);
    Vessel* const vessel = it->second.get();
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

void Plugin::EvolveHistories(Instant const& t) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(t);
  // Integration with a constant step.
  NBodySystem<Barycentric>::Trajectories trajectories;
  // NOTE(egg): This may be too large, vessels that are not new and in the
  // physics bubble or dirty will not be added.
  trajectories.reserve(vessels_.size() - unsynchronized_vessels_.size() +
                       celestials_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_history());
  }
  for (auto const& pair : vessels_) {
    Vessel* const vessel = pair.second.get();
    if (vessel->is_synchronized() &&
        !bubble_.contains(vessel) &&
        !is_dirty(vessel)) {
      trajectories.push_back(vessel->mutable_history());
    }
  }
  VLOG(1) << "Starting the evolution of the histories" << '\n'
          << "from : " << HistoryTime();
  n_body_system_->Integrate(history_integrator_,  // integrator
                            t,                    // tmax
                            Δt_,                  // Δt
                            0,                    // sampling_period
                            false,                // tmax_is_exact
                            trajectories);        // trajectories
  CHECK_GE(HistoryTime(), current_time_);
  VLOG(1) << "Evolved the histories" << '\n'
          << "to   : " << HistoryTime();
}

void Plugin::SynchronizeNewVesselsAndCleanDirtyVessels() {
  VLOG(1) << __FUNCTION__;
  NBodySystem<Barycentric>::Trajectories trajectories;
  trajectories.reserve(celestials_.size() + unsynchronized_vessels_.size() +
                       bubble_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (Vessel* const vessel : unsynchronized_vessels_) {
    if (!bubble_.contains(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  for (Vessel* const vessel : dirty_vessels_) {
    if (!bubble_.contains(vessel) && vessel->is_synchronized()) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  if (!bubble_.empty()) {
    trajectories.push_back(bubble_.mutable_centre_of_mass_trajectory());
  }
  VLOG(1) << "Starting the synchronization of the new vessels"
          << (bubble_.empty() ? "" : " and of the bubble");
  n_body_system_->Integrate(prolongation_integrator_,  // integrator
                            HistoryTime(),             // tmax
                            Δt_,                       // Δt
                            0,                         // sampling_period
                            true,                      // tmax_is_exact
                            trajectories);             // trajectories
  if (!bubble_.empty()) {
    SynchronizeBubbleHistories();
  }
  for (Vessel* const vessel : unsynchronized_vessels_) {
    CHECK(!bubble_.contains(vessel));
    vessel->CreateHistoryAndForkProlongation(
        HistoryTime(),
        vessel->prolongation().last().degrees_of_freedom());
    dirty_vessels_.erase(vessel);
  }
  unsynchronized_vessels_.clear();
  for (Vessel* const vessel : dirty_vessels_) {
    CHECK(!bubble_.contains(vessel));
    vessel->mutable_history()->Append(
        HistoryTime(),
        vessel->prolongation().last().degrees_of_freedom());
  }
  dirty_vessels_.clear();
  VLOG(1) << "Synchronized the new vessels"
          << (bubble_.empty() ? "" : " and the bubble");
}

void Plugin::SynchronizeBubbleHistories() {
  VLOG(1) << __FUNCTION__;
  DegreesOfFreedom<Barycentric> const& centre_of_mass =
      bubble_.centre_of_mass_trajectory().last().degrees_of_freedom();
  for (Vessel* vessel : bubble_.vessels()) {
    RelativeDegreesOfFreedom<Barycentric> const& from_centre_of_mass =
        bubble_.from_centre_of_mass(vessel);
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
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(t);
  NBodySystem<Barycentric>::Trajectories trajectories;
  trajectories.reserve(vessels_.size() + celestials_.size() -
                       bubble_.number_of_vessels() + bubble_.size());
  for (auto const& pair : celestials_) {
    std::unique_ptr<Celestial> const& celestial = pair.second;
    trajectories.push_back(celestial->mutable_prolongation());
  }
  for (auto const& pair : vessels_) {
    Vessel* const vessel = pair.second.get();
    if (!bubble_.contains(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
    }
  }
  if (!bubble_.empty()) {
    trajectories.push_back(bubble_.mutable_centre_of_mass_trajectory());
  }
  VLOG(1) << "Evolving prolongations"
          << (bubble_.empty() ? "" : " and bubble") << '\n'
          << "from : " << trajectories.front()->last().time() << '\n'
          << "to   : " << t;
  n_body_system_->Integrate(prolongation_integrator_,  // integrator
                            t,                         // tmax
                            Δt_,                       // Δt
                            0,                         // sampling_period
                            true,                      // tmax_is_exact
                            trajectories);             // trajectories
  if (!bubble_.empty()) {
    DegreesOfFreedom<Barycentric> const& centre_of_mass =
        bubble_.centre_of_mass_trajectory().last().degrees_of_freedom();
    for (Vessel* vessel : bubble_.vessels()) {
      RelativeDegreesOfFreedom<Barycentric> const& from_centre_of_mass =
          bubble_.from_centre_of_mass(vessel);
      vessel->mutable_prolongation()->Append(
          t,
          centre_of_mass + from_centre_of_mass);
    }
  }
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
    RelativeDegreesOfFreedom<AliceSun> const& from_parent) {
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
  LOG(INFO) << "Initial |{orbit.pos, orbit.vel}| for celestial at index "
            << celestial_index << ": " << from_parent;
  auto const relative =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent));
  LOG(INFO) << "In barycentric coordinates: " << relative;
  Celestial* const celestial = inserted.first->second.get();
  celestial->set_parent(&parent);
  celestial->CreateHistoryAndForkProlongation(
      current_time_,
      parent.history().last().degrees_of_freedom() + relative);
}

void Plugin::EndInitialization() {
  initializing.Flop();
}

void Plugin::UpdateCelestialHierarchy(Index const celestial_index,
                                      Index const parent_index) const {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(celestial_index) << '\n' << NAMED(parent_index);
  CHECK(!initializing);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  auto const it_parent = celestials_.find(parent_index);
  CHECK(it_parent != celestials_.end()) << "No body at index " << parent_index;
  it->second->set_parent(it_parent->second.get());
}

bool Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                Index const parent_index) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(vessel_guid) << '\n' << NAMED(parent_index);
  CHECK(!initializing);
  auto const it = celestials_.find(parent_index);
  CHECK(it != celestials_.end()) << "No body at index " << parent_index;
  Celestial const& parent = *it->second;
  auto inserted = vessels_.emplace(vessel_guid,
                                   std::make_unique<Vessel>(&parent));
  Vessel* const vessel = inserted.first->second.get();
  kept_vessels_.emplace(vessel);
  vessel->set_parent(&parent);
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
  CHECK(!initializing);
  std::unique_ptr<Vessel> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(!vessel->is_initialized())
      << "Vessel with GUID " << vessel_guid << " already has a trajectory";
  LOG(INFO) << "Initial |{orbit.pos, orbit.vel}| for vessel with GUID "
            << vessel_guid << ": " << from_parent;
  RelativeDegreesOfFreedom<Barycentric> const relative =
      PlanetariumRotation().Inverse()(
          kSunLookingGlass.Inverse()(from_parent));
  LOG(INFO) << "In barycentric coordinates: " << relative;
  vessel->CreateProlongation(
      current_time_,
      vessel->parent().history().last().degrees_of_freedom() + relative);
  auto const inserted = unsynchronized_vessels_.emplace(vessel.get());
  CHECK(inserted.second);
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(t) << '\n' << NAMED(planetarium_rotation);
  CHECK(!initializing);
  CHECK_GT(t, current_time_);
  CleanUpVessels();
  bubble_.Prepare(PlanetariumRotation(), current_time_, t);
  if (HistoryTime() + Δt_ < t) {
    // The histories are far enough behind that we can advance them at least one
    // step and reset the prolongations.
    EvolveHistories(t);
    // TODO(egg): I think |!bubble_.empty()| => |has_dirty_vessels()|.
    if (has_unsynchronized_vessels() ||
        has_dirty_vessels() ||
        !bubble_.empty()) {
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
}

RelativeDegreesOfFreedom<AliceSun> Plugin::VesselFromParent(
    GUID const& vessel_guid) const {
  CHECK(!initializing);
  std::unique_ptr<Vessel> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << vessel_guid
                                  << " was not given an initial state";
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      vessel->prolongation().last().degrees_of_freedom() -
      vessel->parent().prolongation().last().degrees_of_freedom();
  RelativeDegreesOfFreedom<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Vessel with GUID " << vessel_guid
          << " is at parent degrees of freedom + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

RelativeDegreesOfFreedom<AliceSun> Plugin::CelestialFromParent(
    Index const celestial_index) const {
  CHECK(!initializing);
  auto const it = celestials_.find(celestial_index);
  CHECK(it != celestials_.end()) << "No body at index " << celestial_index;
  Celestial const& celestial = *it->second;
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      celestial.prolongation().last().degrees_of_freedom() -
      celestial.parent().prolongation().last().degrees_of_freedom();
  RelativeDegreesOfFreedom<AliceSun> const result =
      kSunLookingGlass(PlanetariumRotation()(barycentric_result));
  VLOG(1) << "Celestial at index " << celestial_index
          << " is at parent degrees of freedom + " << barycentric_result
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
          sun_->prolongation().last().degrees_of_freedom().position(),
          sun_world_position,
          Rotation<WorldSun, World>::Identity() * PlanetariumRotation());
  std::unique_ptr<Vessel> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized());
  VLOG(1) << "Rendering a trajectory for the vessel with GUID " << vessel_guid;
  RenderedTrajectory<World> result;
  if (!vessel->is_synchronized()) {
    // TODO(egg): We render neither unsynchronized histories nor prolongations
    // at the moment.
    VLOG(1) << "Returning an empty trajectory";
    return result;
  }
  DegreesOfFreedom<Barycentric> const* initial_state = nullptr;
  DegreesOfFreedom<Barycentric> const* final_state = nullptr;
  std::unique_ptr<Trajectory<Barycentric>> const apparent_trajectory =
      frame.ApparentTrajectory(vessel->history());
  for (Trajectory<Barycentric>::NativeIterator
           it = apparent_trajectory->first();
       !it.at_end();
       ++it) {
    final_state = &it.degrees_of_freedom();
    if (initial_state != nullptr) {
      result.emplace_back(to_world(initial_state->position()),
                          to_world(final_state->position()));
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
  std::unique_ptr<Vessel> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << vessel_guid
                                 << " was not given an initial state";
  auto const to_world =
      AffineMap<Barycentric, World, Length, Rotation>(
          vessel->parent().
              prolongation().last().degrees_of_freedom().position(),
          parent_world_position,
          Rotation<WorldSun, World>::Identity() * PlanetariumRotation());
  return to_world(
      vessel->prolongation().last().degrees_of_freedom().position());
}

Velocity<World> Plugin::VesselWorldVelocity(
      GUID const& vessel_guid,
      Velocity<World> const& parent_world_velocity,
      Time const& parent_rotation_period) const {
  std::unique_ptr<Vessel> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << vessel_guid
                                  << " was not given an initial state";
  Rotation<Barycentric, World> to_world =
      Rotation<WorldSun, World>::Identity() * PlanetariumRotation();
  RelativeDegreesOfFreedom<Barycentric> const relative_to_parent =
      vessel->prolongation().last().degrees_of_freedom() -
      vessel->parent().prolongation().last().degrees_of_freedom();
  AngularVelocity<Barycentric> const world_frame_angular_velocity =
      AngularVelocity<Barycentric>({0 * Radian / Second,
                                    2 * π * Radian / parent_rotation_period,
                                    0 * Radian / Second});
  return to_world(
      (world_frame_angular_velocity *
       relative_to_parent.displacement()) / Radian +
      relative_to_parent.velocity()) + parent_world_velocity;
}
void Plugin::AddVesselToNextPhysicsBubble(
    GUID const& vessel_guid,
    std::vector<IdAndOwnedPart> parts) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel_guid) << '\n' << NAMED(parts);
  std::unique_ptr<Vessel> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  dirty_vessels_.insert(vessel.get());
  bubble_.AddVesselToNext(vessel.get(), std::move(parts));
}

Displacement<World> Plugin::BubbleDisplacementCorrection(
    Position<World> const& sun_world_position) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(sun_world_position);
  VLOG_AND_RETURN(1, bubble_.DisplacementCorrection(PlanetariumRotation(),
                                                    *sun_,
                                                    sun_world_position));
}

Velocity<World> Plugin::BubbleVelocityCorrection(
    Index const reference_body_index) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(reference_body_index);
  auto const found = celestials_.find(reference_body_index);
  CHECK(found != celestials_.end());
  Celestial const& reference_body = *found->second;
  VLOG_AND_RETURN(1, bubble_.VelocityCorrection(PlanetariumRotation(),
                                                reference_body));
}


}  // namespace ksp_plugin
}  // namespace principia
