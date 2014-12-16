#pragma once

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

class PhysicsBubble {
 public:
  using PartCorrespondence = std::pair<Part<World>*, Part<World>*>;

  PhysicsBubble();
  ~PhysicsBubble();

  // Given the vector of common parts produced by |IntrinsicAcceleration|,
  // constructs |*next_physics_bubble_->centre_of_mass_trajectory| and appends
  // degrees of freedom at |current_time_| that conserve the degrees of freedom
  // of the centre of mass of the parts in |common_parts|.
  // |common_parts| must not be null.  |next_physics_bubble_| must not be null.
  // No transfer of ownership.
  void Shift(std::vector<PartCorrespondence> const* const common_parts);

 private:
  using PartIdToOwnedPart = std::map<PartId, std::unique_ptr<Part<World>>>;

  std::map<Vessel* const, std::vector<Part<World>* const>> vessels_;
  PartIdToOwnedPart parts_;
  // TODO(egg): the following six should be |std::optional| when that
  // becomes a thing.
  std::unique_ptr<DegreesOfFreedom<World>> centre_of_mass_;
  std::unique_ptr<Trajectory<Barycentric>> centre_of_mass_trajectory_;
  std::unique_ptr<
      std::map<Vessel const* const,
                Displacement<Barycentric>>> displacements_from_centre_of_mass_;
  std::unique_ptr<
      std::map<Vessel const* const,
                Velocity<Barycentric>>> velocities_from_centre_of_mass_;
  std::unique_ptr<Displacement<World>> displacement_correction_;
  std::unique_ptr<Velocity<World>> velocity_correction_;
  //NOTE(phl): This was the next bubble.
  std::unique_ptr<DegreesOfFreedom<World>> next_centre_of_mass_;
  std::unique_ptr<Trajectory<Barycentric>> next_centre_of_mass_trajectory_;
};

}  // namespace ksp_plugin
}  // namespace principia
