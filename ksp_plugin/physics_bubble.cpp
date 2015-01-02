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

PhysicsBubble::PhysicsBubble()
    : body_() {}

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
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(current_time) << '\n' << NAMED(next_time);
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
      std::vector<PartCorrespondence> const common_parts =
          ComputeCommonParts(*next);
      if (common_parts.empty()) {
        // The current and next set of parts are disjoint, i.e., the next
        // physics bubble is unrelated to the current one.
        RestartNext(current_time, next.get());
      } else {
        Vector<Acceleration, World> const intrinsic_acceleration =
            IntrinsicAcceleration(current_time, next_time, common_parts);
        if (common_parts.size() == next->parts.size() &&
            common_parts.size() == current_->parts.size()) {
          // The set of parts has not changed.
          next->centre_of_mass_trajectory =
              std::move(current_->centre_of_mass_trajectory);
          // TODO(egg): we end up dragging some history along here, we probably
          // should not.
        } else {
          // Parts appeared or were removed from the physics bubble, but the
          // intersection is nonempty.  We fix the degrees of freedom of the
          // centre of mass of the intersection, and we use its measured
          // acceleration as the intrinsic acceleration of the |body_|.
          Shift(planetarium_rotation, current_time, common_parts, next.get());
        }
        // Correct since |World| is currently nonrotating.
        Vector<Acceleration, Barycentric> barycentric_intrinsic_acceleration =
            planetarium_rotation.Inverse()(
                Identity<World, WorldSun>()(intrinsic_acceleration));
        VLOG(1) << NAMED(barycentric_intrinsic_acceleration);
        if (next->centre_of_mass_trajectory->has_intrinsic_acceleration()) {
          next->centre_of_mass_trajectory->clear_intrinsic_acceleration();
        }
        // TODO(egg): this makes the intrinsic acceleration a step function.
        // Might something smoother be better?  We need to be careful not to be
        // one step or half a step in the past though.
        next->centre_of_mass_trajectory->set_intrinsic_acceleration(
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
  CHECK(!empty()) << "Empty bubble";
  if (current_->displacement_correction == nullptr) {
    current_->displacement_correction =
        std::make_unique<Displacement<World>>(
          Identity<WorldSun, World>()(planetarium_rotation(
              current_->centre_of_mass_trajectory->
                  last().degrees_of_freedom().position() -
              reference_celestial.prolongation().
                  last().degrees_of_freedom().position())) +
          reference_celestial_world_position -
              current_->centre_of_mass->position());
  }
  VLOG_AND_RETURN(1, *current_->displacement_correction);
}

Velocity<World> PhysicsBubble::VelocityCorrection(
    PlanetariumRotation const& planetarium_rotation,
    Celestial const& reference_celestial) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(&reference_celestial);
  CHECK(!empty()) << "Empty bubble";
  if (current_->velocity_correction == nullptr) {
    current_->velocity_correction =
        std::make_unique<Velocity<World>>(
            Identity<WorldSun, World>()(planetarium_rotation(
                current_->centre_of_mass_trajectory->
                    last().degrees_of_freedom().velocity() -
                reference_celestial.prolongation().
                    last().degrees_of_freedom().velocity())) -
            current_->centre_of_mass->velocity());
  }
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
  CHECK(!empty()) << "Empty bubble";
  std::vector<Vessel*> vessels;
  for (auto const& pair : current_->vessels) {
    Vessel* const vessel = pair.first;
    vessels.push_back(vessel);
  }
  return vessels;
}

Displacement<Barycentric> const&
PhysicsBubble::displacement_from_centre_of_mass(
    Vessel const* const vessel) const {
  CHECK(!empty()) << "Empty bubble";
  CHECK(current_->displacements_from_centre_of_mass != nullptr);
  auto const it = current_->displacements_from_centre_of_mass->find(vessel);
  CHECK(it != current_->displacements_from_centre_of_mass->end());
  return it->second;
}

Velocity<Barycentric> const&
PhysicsBubble::velocity_from_centre_of_mass(
    Vessel const* const vessel) const {
  CHECK(!empty()) << "Empty bubble";
  CHECK(current_->velocities_from_centre_of_mass != nullptr);
  auto const it = current_->velocities_from_centre_of_mass->find(vessel);
  CHECK(it != current_->velocities_from_centre_of_mass->end());
  return it->second;
}

Trajectory<Barycentric> const&
PhysicsBubble::centre_of_mass_trajectory() const {
  CHECK(!empty()) << "Empty bubble";
  return *current_->centre_of_mass_trajectory;
}

Trajectory<Barycentric>*
PhysicsBubble::mutable_centre_of_mass_trajectory() const {
  CHECK(!empty()) << "Empty bubble";
  return current_->centre_of_mass_trajectory.get();
}

PhysicsBubble::PreliminaryState::PreliminaryState() {}

PhysicsBubble::FullState::FullState(
    PreliminaryState&& preliminary_state)  // NOLINT(build/c++11)
    : PreliminaryState() {
  parts = std::move(preliminary_state.parts);
  vessels = std::move(preliminary_state.vessels);
}

void PhysicsBubble::ComputeNextCentreOfMassWorldDegreesOfFreedom(
    FullState* next) {
  VLOG(1) << __FUNCTION__;
  CHECK_NOTNULL(next);
  DegreesOfFreedom<World>::BarycentreCalculator<Mass> centre_of_mass_calculator;
  for (auto const& id_part : next->parts) {
    std::unique_ptr<Part<World>> const& part = id_part.second;
    centre_of_mass_calculator.Add(part->degrees_of_freedom, part->mass);
  }
  next->centre_of_mass = std::make_unique<DegreesOfFreedom<World>>(
                             centre_of_mass_calculator.Get());
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
  VLOG(1) << NAMED(next->vessels.size());
  for (auto const& vessel_parts : next->vessels) {
    Vessel const* const vessel = vessel_parts.first;
    std::vector<Part<World>* const> const& parts = vessel_parts.second;
    VLOG(1) << NAMED(vessel) << ", " << NAMED(parts.size());
    DegreesOfFreedom<World>::BarycentreCalculator<Mass> vessel_calculator;
    for (auto const part : parts) {
      vessel_calculator.Add(part->degrees_of_freedom, part->mass);
    }
    DegreesOfFreedom<World> const vessel_degrees_of_freedom =
        vessel_calculator.Get();
    Displacement<Barycentric> const displacement_from_centre_of_mass =
        planetarium_rotation.Inverse()(
            Identity<World, WorldSun>()(
                vessel_degrees_of_freedom.position() -
                next->centre_of_mass->position()));
    Velocity<Barycentric> const velocity_from_centre_of_mass =
        planetarium_rotation.Inverse()(
            Identity<World, WorldSun>()(
                vessel_degrees_of_freedom.velocity() -
                next->centre_of_mass->velocity()));
    VLOG(1) << NAMED(displacement_from_centre_of_mass) << ", "
            << NAMED(velocity_from_centre_of_mass);
    next->displacements_from_centre_of_mass->emplace(
        vessel,
        displacement_from_centre_of_mass);
    next->velocities_from_centre_of_mass->emplace(
        vessel,
        velocity_from_centre_of_mass);
  }
}

void PhysicsBubble::RestartNext(Instant const& current_time,
                                FullState* next) {
  VLOG(1) << __FUNCTION__<< '\n' << NAMED(current_time);
  CHECK_NOTNULL(next);
  DegreesOfFreedom<Barycentric>::BarycentreCalculator<Mass> bubble_calculator;
  for (auto const& vessel_parts : next->vessels) {
    Vessel const* vessel = vessel_parts.first;
    std::vector<Part<World>* const> const& parts = vessel_parts.second;
    for (Part<World> const* const part : parts) {
      bubble_calculator.Add(vessel->prolongation().last().degrees_of_freedom(),
                            part->mass);
    }
  }
  next->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(body_);
  next->centre_of_mass_trajectory->Append(current_time,
                                          bubble_calculator.Get());
}

std::vector<PhysicsBubble::PartCorrespondence>
PhysicsBubble::ComputeCommonParts(FullState const& next) {
  VLOG(1) << __FUNCTION__;
  std::vector<PartCorrespondence> common_parts;
  // Most of the time no parts explode.  We reserve accordingly.
  common_parts.reserve(current_->parts.size());
  for (auto it_in_current_parts = current_->parts.cbegin(),
            it_in_next_parts = next.parts.cbegin();
       it_in_current_parts != current_->parts.end() &&
       it_in_next_parts != next.parts.end();) {
    PartId const current_part_id = it_in_current_parts->first;
    PartId const next_part_id = it_in_next_parts->first;
    if (current_part_id < next_part_id) {
      ++it_in_current_parts;
    } else if (next_part_id < current_part_id) {
      ++it_in_next_parts;
    } else {
      std::unique_ptr<Part<World>> const& current_part =
          it_in_current_parts->second;
      std::unique_ptr<Part<World>> const& next_part = it_in_next_parts->second;
      common_parts.emplace_back(current_part.get(), next_part.get());
      ++it_in_current_parts;
      ++it_in_next_parts;
    }
  }
  VLOG_AND_RETURN(1, common_parts);
}

Vector<Acceleration, World> PhysicsBubble::IntrinsicAcceleration(
    Instant const& current_time,
    Instant const& next_time,
    std::vector<PartCorrespondence> const& common_parts) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(current_time) << '\n'
          << NAMED(next_time) << '\n' << NAMED(common_parts);
  CHECK(!common_parts.empty());
  CHECK(current_->velocity_correction != nullptr);
  BarycentreCalculator<Vector<Acceleration, World>, Mass>
      acceleration_calculator;
  Time const δt = next_time - current_time;
  for (auto const& current_next : common_parts) {
    Part<World>* const current_part = current_next.first;
    Part<World>* const next_part = current_next.second;
    acceleration_calculator.Add(
        (next_part->degrees_of_freedom.velocity() -
            (current_part->degrees_of_freedom.velocity() +
             *current_->velocity_correction)) / δt -
        current_part->gravitational_acceleration_to_be_applied_by_ksp,
        // TODO(egg): not sure what we actually want to do here.
        (next_part->mass + current_part->mass) / 2.0);
  }
  VLOG_AND_RETURN(1, acceleration_calculator.Get());
}

void PhysicsBubble::Shift(PlanetariumRotation const& planetarium_rotation,
                          Instant const& current_time,
                          std::vector<PartCorrespondence> const& common_parts,
                          FullState* next) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(current_time) << '\n' << NAMED(common_parts);
  CHECK_NOTNULL(next);
  DegreesOfFreedom<World>::BarycentreCalculator<Mass> current_common_calculator;
  DegreesOfFreedom<World>::BarycentreCalculator<Mass> next_common_calculator;
  for (auto const& current_next : common_parts) {
    Part<World>* const current_part = current_next.first;
    Part<World>* const next_part = current_next.second;
    current_common_calculator.Add(current_part->degrees_of_freedom,
                                  current_part->mass);
    next_common_calculator.Add(next_part->degrees_of_freedom,
                               next_part->mass);
  }
  DegreesOfFreedom<World> const current_common_centre_of_mass =
      current_common_calculator.Get();
  DegreesOfFreedom<World> const next_common_centre_of_mass =
      next_common_calculator.Get();
  // The change in the position of the overall centre of mass resulting from
  // fixing the centre of mass of the intersection.
  auto const a = *next->centre_of_mass - next_common_centre_of_mass;
  auto const c = Identity<World, WorldSun>()(a);
  DegreesOfFreedom<WorldSun> const d = Identity<World, WorldSun>()(next_common_centre_of_mass);

  Displacement<World> const position_change =
      (next->centre_of_mass->position() -
           next_common_centre_of_mass.position()) -
      (current_->centre_of_mass->position() -
           current_common_centre_of_mass.position());
  // The change in the velocity of the overall centre of mass resulting from
  // fixing the velocity of the centre of mass of the intersection.
  Velocity<World> const velocity_change =
      (next->centre_of_mass->velocity() -
           next_common_centre_of_mass.velocity()) -
      (current_->centre_of_mass->velocity() -
           current_common_centre_of_mass.velocity());
  DegreesOfFreedom<Barycentric> const& current_centre_of_mass =
      current_->centre_of_mass_trajectory->last().degrees_of_freedom();
  next->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(body_);
  // Using the identity as the map |World| -> |WorldSun| is valid for
  // velocities too since we assume |World| is currently nonrotating, i.e.,
  // it is stationary with respect to |WorldSun|.
  next->centre_of_mass_trajectory->Append(
      current_time,
      {current_centre_of_mass.position() +
           planetarium_rotation.Inverse()(
               Identity<World, WorldSun>()(position_change)),
       current_centre_of_mass.velocity() +
           planetarium_rotation.Inverse()(
               Identity<World, WorldSun>()(velocity_change))});
}

}  // namespace ksp_plugin
}  // namespace principia
