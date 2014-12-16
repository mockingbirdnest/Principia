#pragma once

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "geometry/rotation.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

class PhysicsBubble {
 public:
  using PartCorrespondence = std::pair<Part<World>*, Part<World>*>;
  using PlanetariumRotationXXX = geometry::Rotation<Barycentric, WorldSun>;

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  PhysicsBubble() = default;
  ~PhysicsBubble() = default;

  // If |next_| is not null, computes the world centre of mass, trajectory
  // (including intrinsic acceleration) of |*next_|. Moves |next_| into
  // |current_|.
  // TODO(phl): Document the parameters!
  void Prepare(PlanetariumRotationXXX const& planetarium_rotation,
               Instant const& current_time,
               Instant const& next_time);

 private:
  // Computes the world degrees of freedom of the centre of mass of
  // |next_| using the contents of |next_->parts|.  |next_| must not be null.
  void ComputeNextCentreOfMassWorldDegreesOfFreedom();

  // Computes |next_->displacements_from_centre_of_mass| and
  // |next_->velocities_from_centre_of_mass|.  |next_| must not be null.
  void ComputeNextVesselOffsets(
      PlanetariumRotationXXX const& planetarium_rotation);

  // Creates |next_->centre_of_mass_trajectory| and appends to it the barycentre
  // of the degrees of freedom of the vessels in |next_->vessels|.  There is no
  // intrinsic acceleration.  |next_| must not be null.
  void RestartNext(Instant const& current_time);

  // Returns the intrinsic acceleration measured on the parts that are common to
  // the current and next bubbles.  Stores a pair of pointers to parts
  // (current, next) in |common_parts| for all parts common to the current and
  // next bubbles.  |common_parts| must not be null.  |*common_parts| must be
  // empty.  No transfer of ownership.
  Vector<Acceleration, World> IntrinsicAcceleration(
      Instant const& current_time,
      Instant const& next_time,
      std::vector<PartCorrespondence>* const common_parts);

  // Given the vector of common parts produced by |IntrinsicAcceleration|,
  // constructs |*next_->centre_of_mass_trajectory| and appends degrees of
  // freedom at |current_time| that conserve the degrees of freedom of the
  // centre of mass of the parts in |common_parts|. |common_parts| must not be
  // null.  |next_| must not be null.  No transfer of ownership.
  void Shift(PlanetariumRotationXXX const& planetarium_rotation,
             Instant const& current_time,
             std::vector<PartCorrespondence> const* const common_parts);

  using PartIdToOwnedPart = std::map<PartId, std::unique_ptr<Part<World>>>;

  struct State {
    std::map<Vessel* const, std::vector<Part<World>* const>> vessels;
    PartIdToOwnedPart parts;
    // TODO(egg): the following six should be |std::optional| when that
    // becomes a thing.
    std::unique_ptr<DegreesOfFreedom<World>> centre_of_mass;
    std::unique_ptr<Trajectory<Barycentric>> centre_of_mass_trajectory;
    std::unique_ptr<
        std::map<Vessel const* const,
                 Displacement<Barycentric>>> displacements_from_centre_of_mass;
    std::unique_ptr<
        std::map<Vessel const* const,
                 Velocity<Barycentric>>> velocities_from_centre_of_mass;
    std::unique_ptr<Displacement<World>> displacement_correction;  // Only for current_.
    std::unique_ptr<Velocity<World>> velocity_correction;  // Only for current_.
  };

  //TODO(phl): See what's used exactly for next_.
  std::unique_ptr<State> current_;
  std::unique_ptr<State> next_;

  MasslessBody const bubble_body_;
};

}  // namespace ksp_plugin
}  // namespace principia
