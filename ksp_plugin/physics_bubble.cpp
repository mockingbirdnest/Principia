#include "ksp_plugin/physics_bubble.hpp"

#include <vector>

#include "geometry/identity.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"

using principia::geometry::Identity;

namespace principia {
namespace ksp_plugin {

void PhysicsBubble::Shift(
    std::vector<PartCorrespondence> const* const common_parts) {
  VLOG(1) << __FUNCTION__ << '\n'<< NAMED(common_parts);
  CHECK_NOTNULL(common_parts);
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
      (next_centre_of_mass_->position -
           next_common_centre_of_mass.position) -
      (centre_of_mass_->position -
           current_common_centre_of_mass.position);
  // The change in the velocity of the overall centre of mass resulting from
  // fixing the velocity of the centre of mass of the intersection.
  Velocity<World> const velocity_change =
      (next_centre_of_mass_->velocity -
           next_common_centre_of_mass.velocity) -
      (centre_of_mass_->velocity -
           current_common_centre_of_mass.velocity);
  DegreesOfFreedom<Barycentric> const& current_centre_of_mass =
      centre_of_mass_trajectory_->last().degrees_of_freedom();
  next_centre_of_mass_trajectory_ =
      std::make_unique<Trajectory<Barycentric>>(bubble_body_);
  // Using the identity as the map |World| -> |WorldSun| is valid for
  // velocities too since we assume |World| is currently nonrotating, i.e.,
  // it is stationary with respect to |WorldSun|.
  next_centre_of_mass_trajectory_->Append(
      current_time_,
      {current_centre_of_mass.position +
           PlanetariumRotation().Inverse()(
               Identity<World, WorldSun>()(position_change)),
       current_centre_of_mass.velocity +
           PlanetariumRotation().Inverse()(
               Identity<World, WorldSun>()(velocity_change))});
}


}  // namespace ksp_plugin
}  // namespace principia
