#include "ksp_plugin/physics_bubble.hpp"

#include <map>
#include <utility>
#include <vector>

#include "base/macros.hpp"
#include "base/unique_ptr_logging.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/identity.hpp"
#include "glog/stl_logging.h"
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::BarycentreCalculator;
using principia::geometry::Identity;
using principia::quantities::Time;

namespace principia {
namespace ksp_plugin {

void PhysicsBubble::AddVesselToNext(Vessel* vessel,
                                    std::vector<IdAndOwnedPart> parts) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel) << '\n' << NAMED(parts);
  if (next_ == nullptr) {
    next_ = std::make_unique<PreliminaryState>();
  }
  auto const inserted_vessel =
      next_->vessels.emplace(vessel, std::vector<Part<World>* const>());
  CHECK(inserted_vessel.second);
  std::vector<Part<World>* const>* const vessel_parts =
      &inserted_vessel.first->second;
  for (IdAndOwnedPart& id_part : parts) {
    PartId const id = id_part.first;
    std::unique_ptr<Part<World>> const& part = id_part.second;
    VLOG(1) << "Inserting {id, part}" << '\n' << NAMED(id) << '\n'
            << NAMED(*part);
    auto const inserted_part =
        next_->parts.insert(std::move(id_part));
    CHECK(inserted_part.second) << id;
    VLOG(1) << "Part is at: " << inserted_part.first->second;
    vessel_parts->push_back(inserted_part.first->second.get());
  }
}

void PhysicsBubble::Prepare(PlanetariumRotation const& planetarium_rotation,
                            Instant const& current_time,
                            Instant const& next_time) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(next_time);
  std::unique_ptr<FullState> next;
  if (next_ != nullptr) {
    next = std::make_unique<FullState>(std::move(*next_));
    next_.reset();
    ComputeNextCentreOfMassWorldDegreesOfFreedom(next.get());
    ComputeNextVesselOffsets(planetarium_rotation, next.get());
    if (current_ == nullptr) {
      // There was no physics bubble.
      RestartNext(current_time, next.get());
    } else {
      // The IDs of the parts that are both in the current and in the next
      // physics bubble.
      auto const common_parts =
          std::make_unique<
              std::vector<std::pair<Part<World>*, Part<World>*>>>();
      Vector<Acceleration, World> const intrinsic_acceleration =
          IntrinsicAcceleration(current_time, next_time, common_parts.get());
      if (common_parts->empty()) {
        // The current and next set of parts are disjoint, i.e., the next
        // physics bubble is unrelated to the current one.
        RestartNext(current_time, next.get());
      } else {
        if (common_parts->size() == next_->parts.size() &&
            common_parts->size() == current_->parts.size()) {
          // The set of parts has not changed.
          next->centre_of_mass_trajectory =
              std::move(current_->centre_of_mass_trajectory);
          // TODO(egg): we end up dragging some history along here, we probably
          // should not.
        } else {
          // Parts appeared or were removed from the physics bubble, but the
          // intersection is nonempty.  We fix the degrees of freedom of the
          // centre of mass of the intersection, and we use its measured
          // acceleration as the intrinsic acceleration of the |bubble_body_|.
          Shift(planetarium_rotation,
                current_time,
                common_parts.get(),
                next.get());
        }
        // Correct since |World| is currently nonrotating.
        Vector<Acceleration, Barycentric> barycentric_intrinsic_acceleration =
            planetarium_rotation.Inverse()(
                Identity<World, WorldSun>()(intrinsic_acceleration));
        VLOG(1) << NAMED(barycentric_intrinsic_acceleration);
        if (next->centre_of_mass_trajectory->
                has_intrinsic_acceleration()) {
          next->centre_of_mass_trajectory->
              clear_intrinsic_acceleration();
        }
        // TODO(egg): this makes the intrinsic acceleration a step function.
        // Might something smoother be better?  We need to be careful not to be
        // one step or half a step in the past though.
        next->centre_of_mass_trajectory->
              set_intrinsic_acceleration(
                  [barycentric_intrinsic_acceleration](Instant const& t) {
                    return barycentric_intrinsic_acceleration;
                  });
      }
    }
  }
  current_ = std::move(next);
  CHECK(next_ == nullptr);
  VLOG_IF(1, current_ == nullptr) << "No physics bubble";
  VLOG_IF(1, current_ != nullptr)
      << "Bubble will be integrated from: "
      << current_->centre_of_mass_trajectory->last().degrees_of_freedom();
}

Displacement<World> PhysicsBubble::DisplacementCorrection(
    PlanetariumRotation const& planetarium_rotation,
    Celestial const& reference_celestial,
    Position<World> const& reference_celestial_world_position) const {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(&reference_celestial) << '\n'
          << NAMED(reference_celestial_world_position);
  CHECK(!empty());
  CHECK(current_->displacement_correction == nullptr);
  current_->displacement_correction =
      std::make_unique<Displacement<World>>(
        Identity<WorldSun, World>()(planetarium_rotation(
            current_->centre_of_mass_trajectory->
                last().degrees_of_freedom().position -
            reference_celestial.prolongation().
                last().degrees_of_freedom().position)) +
        reference_celestial_world_position -
            current_->centre_of_mass->position);
  VLOG_AND_RETURN(1, *current_->displacement_correction);
}

Velocity<World> PhysicsBubble::VelocityCorrection(
    PlanetariumRotation const& planetarium_rotation,
    Celestial const& reference_celestial) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(&reference_celestial);
  CHECK(!empty());
  CHECK(current_->velocity_correction == nullptr);
  current_->velocity_correction =
      std::make_unique<Velocity<World>>(
          Identity<WorldSun, World>()(planetarium_rotation(
              current_->centre_of_mass_trajectory->
                  last().degrees_of_freedom().velocity -
              reference_celestial.prolongation().
                  last().degrees_of_freedom().velocity)) -
          current_->centre_of_mass->velocity);
  VLOG_AND_RETURN(1, *current_->velocity_correction);
}

bool PhysicsBubble::empty() const {
  return current_ == nullptr;
}

std::size_t PhysicsBubble::size() const {
  return empty() ? 0 : 1;
}

std::size_t PhysicsBubble::number_of_vessels() const {
  if (empty()) {
    return 0;
  } else {
    return current_->vessels.size();
  }
}

bool PhysicsBubble::contains(Vessel* const vessel) const {
  return !empty() &&
         current_->vessels.find(vessel) != current_->vessels.end();
}

std::vector<Vessel*> PhysicsBubble::vessels() const {
  CHECK(!empty());
  std::vector<Vessel*> vessels;
  for (auto const& pair : current_->vessels) {
    Vessel* const vessel = pair.first;
    vessels.push_back(vessel);
  }
  return vessels;
}

Displacement<Barycentric> const&
PhysicsBubble::displacements_from_centre_of_mass(
    Vessel const* const vessel) const {
  CHECK(!empty());
  CHECK(current_->displacements_from_centre_of_mass != nullptr);
  auto const it = current_->displacements_from_centre_of_mass->find(vessel);
  CHECK(it != current_->displacements_from_centre_of_mass->end());
  return it->second;
}

Velocity<Barycentric> const&
PhysicsBubble::velocities_from_centre_of_mass(
    Vessel const* const vessel) const {
  CHECK(!empty());
  CHECK(current_->velocities_from_centre_of_mass != nullptr);
  auto const it = current_->velocities_from_centre_of_mass->find(vessel);
  CHECK(it != current_->velocities_from_centre_of_mass->end());
  return it->second;
}

Trajectory<Barycentric> const&
PhysicsBubble::centre_of_mass_trajectory() const {
  CHECK(!empty());
  return *current_->centre_of_mass_trajectory;
}

Trajectory<Barycentric>*
PhysicsBubble::mutable_centre_of_mass_trajectory() const {
  CHECK(!empty());
  return current_->centre_of_mass_trajectory.get();
}

PhysicsBubble::PreliminaryState::PreliminaryState() {}

PhysicsBubble::FullState::FullState(PreliminaryState&& preliminary_state)
    : PreliminaryState() {
  parts = std::move(preliminary_state.parts);
  vessels = std::move(preliminary_state.vessels);
}

void PhysicsBubble::ComputeNextCentreOfMassWorldDegreesOfFreedom(
    FullState* next) {
  VLOG(1) << __FUNCTION__;
  CHECK_NOTNULL(next);
  std::vector<DegreesOfFreedom<World>> part_degrees_of_freedom;
  part_degrees_of_freedom.reserve(next_->parts.size());
  std::vector<Mass> part_masses;
  part_masses.reserve(next_->parts.size());
  for (auto const& id_part : next_->parts) {
    std::unique_ptr<Part<World>> const& part = id_part.second;
    part_degrees_of_freedom.push_back(part->degrees_of_freedom);
    part_masses.push_back(part->mass);
  }
  next->centre_of_mass =
      std::make_unique<DegreesOfFreedom<World>>(
          physics::Barycentre(part_degrees_of_freedom, part_masses));
  VLOG(1) << NAMED(*next->centre_of_mass);
}

void PhysicsBubble::ComputeNextVesselOffsets(
    PlanetariumRotation const& planetarium_rotation,
    FullState* next) {
  VLOG(1) << __FUNCTION__;
  CHECK_NOTNULL(next);
  next->displacements_from_centre_of_mass =
      std::make_unique<std::map<Vessel const* const,
                                Displacement<Barycentric>>>();
  next->velocities_from_centre_of_mass =
      std::make_unique<std::map<Vessel const* const,
                                Velocity<Barycentric>>>();
  VLOG(1) << NAMED(next_->vessels.size());
  for (auto const& vessel_parts : next_->vessels) {
    Vessel const* const vessel = vessel_parts.first;
    std::vector<Part<World>* const> const& parts = vessel_parts.second;
    VLOG(1) << NAMED(vessel) << ", " << NAMED(parts.size());
    std::vector<DegreesOfFreedom<World>> part_degrees_of_freedom;
    std::vector<Mass> part_masses;
    part_degrees_of_freedom.reserve(parts.size());
    part_masses.reserve(parts.size());
    for (auto const part : parts) {
      part_degrees_of_freedom.emplace_back(part->degrees_of_freedom);
      part_masses.emplace_back(part->mass);
    }
    DegreesOfFreedom<World> const vessel_degrees_of_freedom =
        physics::Barycentre(part_degrees_of_freedom, part_masses);
    Displacement<Barycentric> const displacement_from_centre_of_mass =
        planetarium_rotation.Inverse()(
            Identity<World, WorldSun>()(
                vessel_degrees_of_freedom.position -
                next->centre_of_mass->position));
    Velocity<Barycentric> const velocity_from_centre_of_mass =
        planetarium_rotation.Inverse()(
            Identity<World, WorldSun>()(
                vessel_degrees_of_freedom.velocity -
                next->centre_of_mass->velocity));
    VLOG(1) << NAMED(displacement_from_centre_of_mass) << ", "
            << NAMED(velocity_from_centre_of_mass);
    next->displacements_from_centre_of_mass->emplace(
        vessel,
        displacement_from_centre_of_mass);
    next->velocities_from_centre_of_mass->emplace(
        vessel,
        velocity_from_centre_of_mass);
  }}

void PhysicsBubble::RestartNext(Instant const& current_time,
                                FullState* next) {
  VLOG(1) << __FUNCTION__;
  CHECK_NOTNULL(next);
  std::vector<DegreesOfFreedom<Barycentric>> vessel_degrees_of_freedom;
  vessel_degrees_of_freedom.reserve(next_->vessels.size());
  std::vector<Mass> vessel_masses;
  vessel_masses.reserve(next_->vessels.size());
  for (auto const& vessel_parts : next_->vessels) {
    Vessel const* vessel = vessel_parts.first;
    std::vector<Part<World>* const> const& parts = vessel_parts.second;
    vessel_degrees_of_freedom.push_back(
        vessel->prolongation().last().degrees_of_freedom());
    vessel_masses.push_back(Mass());
    for (Part<World> const* const part : parts) {
      vessel_masses.back() += part->mass;
    }
  }
  next->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(bubble_body_);
  next->centre_of_mass_trajectory->Append(
      current_time,
      physics::Barycentre(vessel_degrees_of_freedom, vessel_masses));
}

Vector<Acceleration, World> PhysicsBubble::IntrinsicAcceleration(
    Instant const& current_time,
    Instant const& next_time,
    std::vector<PartCorrespondence>* const common_parts) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(next_time) << '\n' << NAMED(common_parts);
  CHECK_NOTNULL(common_parts);
  CHECK(common_parts->empty());
  CHECK(current_->velocity_correction != nullptr);
  // Most of the time no parts explode.  We reserve accordingly.
  common_parts->reserve(current_->parts.size());
  BarycentreCalculator<Vector<Acceleration, World>, Mass>
      acceleration_calculator;
  Time const δt = next_time - current_time;
  for (auto it_in_current_parts = current_->parts.cbegin(),
            it_in_next_parts = next_->parts.cbegin();
       it_in_current_parts != current_->parts.end() &&
       it_in_next_parts != next_->parts.end();) {
    PartId current_part_id = it_in_current_parts->first;
    PartId next_part_id = it_in_next_parts->first;
    if (current_part_id < next_part_id) {
      ++it_in_current_parts;
    } else if (next_part_id < current_part_id) {
      ++it_in_next_parts;
    } else {
      std::unique_ptr<Part<World>> const& current_part =
          it_in_current_parts->second;
      std::unique_ptr<Part<World>> const& next_part = it_in_next_parts->second;
      common_parts->emplace_back(current_part.get(), next_part.get());
      acceleration_calculator.Add(
          (next_part->degrees_of_freedom.velocity -
           (current_part->degrees_of_freedom.velocity +
            *current_->velocity_correction)) / δt -
          current_part->gravitational_acceleration_to_be_applied_by_ksp,
          // TODO(egg): not sure what we actually want to do here.
          (next_part->mass + current_part->mass) / 2.0);
      ++it_in_current_parts;
      ++it_in_next_parts;
    }
  }
  VLOG(1) << NAMED(*common_parts);
  VLOG_AND_RETURN(1, acceleration_calculator.Get());
}

void PhysicsBubble::Shift(
    PlanetariumRotation const& planetarium_rotation,
    Instant const& current_time,
    std::vector<PartCorrespondence> const* const common_parts,
    FullState* next) {
  VLOG(1) << __FUNCTION__ << '\n'<< NAMED(common_parts);
  CHECK_NOTNULL(common_parts);
  CHECK_NOTNULL(next);
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
      (next->centre_of_mass->position -
           next_common_centre_of_mass.position) -
      (current_->centre_of_mass->position -
           current_common_centre_of_mass.position);
  // The change in the velocity of the overall centre of mass resulting from
  // fixing the velocity of the centre of mass of the intersection.
  Velocity<World> const velocity_change =
      (next->centre_of_mass->velocity -
           next_common_centre_of_mass.velocity) -
      (current_->centre_of_mass->velocity -
           current_common_centre_of_mass.velocity);
  DegreesOfFreedom<Barycentric> const& current_centre_of_mass =
      current_->
          centre_of_mass_trajectory->last().degrees_of_freedom();
  next->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(bubble_body_);
  // Using the identity as the map |World| -> |WorldSun| is valid for
  // velocities too since we assume |World| is currently nonrotating, i.e.,
  // it is stationary with respect to |WorldSun|.
  next->centre_of_mass_trajectory->Append(
      current_time,
      {current_centre_of_mass.position +
           planetarium_rotation.Inverse()(
               Identity<World, WorldSun>()(position_change)),
       current_centre_of_mass.velocity +
           planetarium_rotation.Inverse()(
               Identity<World, WorldSun>()(velocity_change))});
}

}  // namespace ksp_plugin
}  // namespace principia
