#include "ksp_plugin/physics_bubble.hpp"

#include <vector>

#include "base/macros.hpp"
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

void PhysicsBubble::Prepare(PlanetariumRotationXXX const& planetarium_rotation,
                            Instant const& next_time) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(next_time);
  if (next_ != nullptr) {
    ComputeNextCentreOfMassWorldDegreesOfFreedom();
    ComputeNextVesselOffsets(planetarium_rotation);
    if (current_ == nullptr) {
      // There was no physics bubble.
      RestartNext();
    } else {
      // The IDs of the parts that are both in the current and in the next
      // physics bubble.
      auto const common_parts =
          std::make_unique<
              std::vector<std::pair<Part<World>*, Part<World>*>>>();
      Vector<Acceleration, World> const intrinsic_acceleration =
          IntrinsicAcceleration(next_time, common_parts.get());
      if (common_parts->empty()) {
        // The current and next set of parts are disjoint, i.e., the next
        // physics bubble is unrelated to the current one.
        RestartNext();
      } else {
        if (common_parts->size() == next_->parts.size() &&
            common_parts->size() == current_->parts.size()) {
          // The set of parts has not changed.
          next_->centre_of_mass_trajectory =
              std::move(current_->centre_of_mass_trajectory);
          // TODO(egg): we end up dragging some history along here, we probably
          // should not.
        } else {
          // Parts appeared or were removed from the physics bubble, but the
          // intersection is nonempty.  We fix the degrees of freedom of the
          // centre of mass of the intersection, and we use its measured
          // acceleration as the intrinsic acceleration of the |bubble_body_|.
          Shift(planetarium_rotation, common_parts.get());
        }
        // Correct since |World| is currently nonrotating.
        Vector<Acceleration, Barycentric> barycentric_intrinsic_acceleration =
            planetarium_rotation.Inverse()(
                Identity<World, WorldSun>()(intrinsic_acceleration));
        VLOG(1) << NAMED(barycentric_intrinsic_acceleration);
        if (next_->centre_of_mass_trajectory->
                has_intrinsic_acceleration()) {
          next_->centre_of_mass_trajectory->
              clear_intrinsic_acceleration();
        }
        // TODO(egg): this makes the intrinsic acceleration a step function.
        // Might something smoother be better?  We need to be careful not to be
        // one step or half a step in the past though.
        next_->centre_of_mass_trajectory->
              set_intrinsic_acceleration(
                  [barycentric_intrinsic_acceleration](Instant const& t) {
                    return barycentric_intrinsic_acceleration;
                  });
      }
    }
  }
  current_ = std::move(next_);
  VLOG_IF(1, current_ == nullptr) << "No physics bubble";
  VLOG_IF(1, current_ != nullptr)
      << "Bubble will be integrated from: "
      << current_->centre_of_mass_trajectory->last().degrees_of_freedom();
}

void PhysicsBubble::ComputeNextCentreOfMassWorldDegreesOfFreedom() {
  VLOG(1) << __FUNCTION__;
  CHECK(next_ != nullptr);
  std::vector<DegreesOfFreedom<World>> part_degrees_of_freedom;
  part_degrees_of_freedom.reserve(next_->parts.size());
  std::vector<Mass> part_masses;
  part_masses.reserve(next_->parts.size());
  for (auto const& id_part : next_->parts) {
    std::unique_ptr<Part<World>> const& part = id_part.second;
    part_degrees_of_freedom.push_back(part->degrees_of_freedom);
    part_masses.push_back(part->mass);
  }
  next_->centre_of_mass =
      std::make_unique<DegreesOfFreedom<World>>(
          physics::Barycentre(part_degrees_of_freedom, part_masses));
  VLOG(1) << NAMED(*next_->centre_of_mass);
}

void PhysicsBubble::ComputeNextVesselOffsets(
    PlanetariumRotationXXX const& planetarium_rotation) {
  VLOG(1) << __FUNCTION__;
  CHECK(next_ != nullptr);
  next_->displacements_from_centre_of_mass =
      std::make_unique<std::map<Vessel const* const,
                                Displacement<Barycentric>>>();
  next_->velocities_from_centre_of_mass =
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
                next_->centre_of_mass->position));
    Velocity<Barycentric> const velocity_from_centre_of_mass =
        planetarium_rotation.Inverse()(
            Identity<World, WorldSun>()(
                vessel_degrees_of_freedom.velocity -
                next_->centre_of_mass->velocity));
    VLOG(1) << NAMED(displacement_from_centre_of_mass) << ", "
            << NAMED(velocity_from_centre_of_mass);
    next_->displacements_from_centre_of_mass->emplace(
        vessel,
        displacement_from_centre_of_mass);
    next_->velocities_from_centre_of_mass->emplace(
        vessel,
        velocity_from_centre_of_mass);
  }}

void PhysicsBubble::RestartNext() {
  VLOG(1) << __FUNCTION__;
  CHECK(next_ != nullptr);
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
  next_->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(bubble_body_);
  next_->centre_of_mass_trajectory->Append(
      current_time_,
      physics::Barycentre(vessel_degrees_of_freedom, vessel_masses));
}

Vector<Acceleration, World> PhysicsBubble::IntrinsicAcceleration(
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
  Time const δt = next_time - current_time_;
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
    PlanetariumRotationXXX const& planetarium_rotation,
    std::vector<PartCorrespondence> const* const common_parts) {
  VLOG(1) << __FUNCTION__ << '\n'<< NAMED(common_parts);
  CHECK_NOTNULL(common_parts);
  CHECK(next_ != nullptr);
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
      (next_->centre_of_mass->position -
           next_common_centre_of_mass.position) -
      (current_->centre_of_mass->position -
           current_common_centre_of_mass.position);
  // The change in the velocity of the overall centre of mass resulting from
  // fixing the velocity of the centre of mass of the intersection.
  Velocity<World> const velocity_change =
      (next_->centre_of_mass->velocity -
           next_common_centre_of_mass.velocity) -
      (current_->centre_of_mass->velocity -
           current_common_centre_of_mass.velocity);
  DegreesOfFreedom<Barycentric> const& current_centre_of_mass =
      current_->
          centre_of_mass_trajectory->last().degrees_of_freedom();
  next_->centre_of_mass_trajectory =
      std::make_unique<Trajectory<Barycentric>>(bubble_body_);
  // Using the identity as the map |World| -> |WorldSun| is valid for
  // velocities too since we assume |World| is currently nonrotating, i.e.,
  // it is stationary with respect to |WorldSun|.
  next_->centre_of_mass_trajectory->Append(
      current_time_,
      {current_centre_of_mass.position +
           planetarium_rotation.Inverse()(
               Identity<World, WorldSun>()(position_change)),
       current_centre_of_mass.velocity +
           planetarium_rotation.Inverse()(
               Identity<World, WorldSun>()(velocity_change))});
}

}  // namespace ksp_plugin
}  // namespace principia
