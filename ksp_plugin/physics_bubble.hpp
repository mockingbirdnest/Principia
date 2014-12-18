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
  //TODO(phl): Move to some common place.
  using PartIdToOwnedPart = std::map<PartId, std::unique_ptr<Part<World>>>;
  using IdAndOwnedPart = PartIdToOwnedPart::value_type;
  using Index = int;
  using PartCorrespondence = std::pair<Part<World>*, Part<World>*>;
  using PlanetariumRotation = geometry::Rotation<Barycentric, WorldSun>;

  PhysicsBubble() = default;
  ~PhysicsBubble() = default;

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  // Creates |next_physics_bubble_| if it is null.  Adds the vessel with GUID
  // |vessel_guid| to |next_physics_bubble_->vessels| with a list of pointers to
  // the |Part|s in |parts|.  Merges |parts| into |next_physics_bubble_->parts|.
  // Adds the vessel to |dirty_vessels_|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept.  The
  // vessel with GUID |vessel_guid| must not already be in
  // |next_physics_bubble_->vessels|.  |parts| must not contain a |PartId|
  // already in |next_physics_bubble_->parts|.
  void PhysicsBubble::AddVesselToNextPhysicsBubble(
      Vessel* vessel,
      std::vector<IdAndOwnedPart> parts);

  // If |next_| is not null, computes the world centre of mass, trajectory
  // (including intrinsic acceleration) of |*next_|. Moves |next_| into
  // |current_|.
  // TODO(phl): Document the parameters!
  void Prepare(PlanetariumRotation const& planetarium_rotation,
               Instant const& current_time,
               Instant const& next_time);

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  Displacement<World> PhysicsBubble::DisplacementCorrection(
      PlanetariumRotation const& planetarium_rotation,
      Celestial const& reference_celestial,
      Position<World> const& reference_celestial_world_position) const;

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  Velocity<World> PhysicsBubble::VelocityCorrection(
      PlanetariumRotation const& planetarium_rotation,
      Celestial const& reference_celestial) const;

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  // Returns |current_physics_bubble_ != nullptr|.
  bool empty() const;

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  // Returns 1 if |has_physics_bubble()|, 0 otherwise.
  std::size_t size() const;

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  // Returns |current_physics_bubble_->vessels.size()|, or 0 if
  // |current_physics_bubble_| is null.
  std::size_t number_of_vessels() const;

  //TODO(phl): Fix \o/ ALL \o/ the comments.
  // Returns true if, and only if, |vessel| is in
  // |current_physics_bubble_->vessels|.  |current_physics_bubble_| may be null,
  // in that case, returns false.
  bool is_in_physics_bubble(Vessel* const vessel) const;

  //TODO(phl): comments.
  std::vector<Vessel const*> vessels() const;
  Displacement<Barycentric> const& displacements_from_centre_of_mass(
      Vessel* const vessel) const;
  Velocity<Barycentric> const& velocities_from_centre_of_mass(
      Vessel* const vessel) const;
  Trajectory<Barycentric> const& centre_of_mass_trajectory() const;
  Trajectory<Barycentric>* mutable_centre_of_mass_trajectory() const;

private:
  struct PreliminaryState {
    std::map<Vessel* const, std::vector<Part<World>* const>> vessels;
    PartIdToOwnedPart parts;
  };

  struct FullState : public PreliminaryState {
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

  // Computes the world degrees of freedom of the centre of mass of
  // |next_| using the contents of |next_->parts|.  |next_| must not be null.
  void ComputeNextCentreOfMassWorldDegreesOfFreedom(FullState* next);

  // Computes |next_->displacements_from_centre_of_mass| and
  // |next_->velocities_from_centre_of_mass|.  |next_| must not be null.
  void ComputeNextVesselOffsets(
      PlanetariumRotation const& planetarium_rotation,
      FullState* next);

  // Creates |next_->centre_of_mass_trajectory| and appends to it the barycentre
  // of the degrees of freedom of the vessels in |next_->vessels|.  There is no
  // intrinsic acceleration.  |next_| must not be null.
  void RestartNext(Instant const& current_time,
                   FullState* next);

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
  void Shift(PlanetariumRotation const& planetarium_rotation,
             Instant const& current_time,
             std::vector<PartCorrespondence> const* const common_parts,
             FullState* next);

  std::unique_ptr<FullState> current_;
  std::unique_ptr<PreliminaryState> next_;

  MasslessBody const bubble_body_;
};

}  // namespace ksp_plugin
}  // namespace principia
