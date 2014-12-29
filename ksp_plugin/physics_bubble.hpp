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
  using PlanetariumRotation = geometry::Rotation<Barycentric, WorldSun>;

  PhysicsBubble();
  ~PhysicsBubble() = default;

  // Creates |next_| if it is null.  Adds the |vessel| to |next_->vessels| with
  // a list of pointers to the Parts in |parts|.  Merges |parts| into
  // |next_->parts|.  The |vessel| must not already be in |next_->vessels|.
  // |parts| must not contain a |PartId| already in |next_->parts|.
  void AddVesselToNext(Vessel* vessel, std::vector<IdAndOwnedPart> parts);

  // If |next_| is not null, computes the world centre of mass, trajectory
  // (including intrinsic acceleration) of |*next_|. Moves |next_| into
  // |current_|.
  // TODO(phl): Document the parameters!
  void Prepare(PlanetariumRotation const& planetarium_rotation,
               Instant const& current_time,
               Instant const& next_time);

  // Computes and returns |current_->displacement_correction|.  This is the
  // |World| shift to be applied to the bubble in order for it to be in the
  // correct position.
  Displacement<World> DisplacementCorrection(
      PlanetariumRotation const& planetarium_rotation,
      Celestial const& reference_celestial,
      Position<World> const& reference_celestial_world_position) const;

  // Computes and returns |current_->velocity_correction|.  This is the |World|
  // shift to be applied to the physics bubble in order for it to have the
  // correct velocity.
  Velocity<World> VelocityCorrection(
      PlanetariumRotation const& planetarium_rotation,
      Celestial const& reference_celestial) const;

  // Returns |current_ != nullptr|.
  bool empty() const;

  // Returns 0 if |empty()|, 1 otherwise.
  std::size_t size() const;

  // Returns |current_->vessels.size()|, or 0 if |empty()|.
  std::size_t number_of_vessels() const;

  // Returns true if, and only if, |vessel| is in |current_->vessels|.
  // |current_| may be null, in that case, returns false.
  bool contains(Vessel* const vessel) const;

  // Selectors for the data in |current_|.
  std::vector<Vessel*> vessels() const;
  Displacement<Barycentric> const& displacement_from_centre_of_mass(
      Vessel const* const vessel) const;
  Velocity<Barycentric> const& velocity_from_centre_of_mass(
      Vessel const* const vessel) const;
  Trajectory<Barycentric> const& centre_of_mass_trajectory() const;
  Trajectory<Barycentric>* mutable_centre_of_mass_trajectory() const;

 private:
  using PartCorrespondence = std::pair<Part<World>*, Part<World>*>;

  struct PreliminaryState {
    PreliminaryState();
    std::map<Vessel* const, std::vector<Part<World>* const>> vessels;
    PartIdToOwnedPart parts;
  };

  struct FullState : public PreliminaryState {
    explicit FullState(
        PreliminaryState&& preliminary_state);  // NOLINT(build/c++11)

    // TODO(egg): these fields should be |std::optional| when that becomes a
    // thing.
    std::unique_ptr<DegreesOfFreedom<World>> centre_of_mass;
    std::unique_ptr<Trajectory<Barycentric>> centre_of_mass_trajectory;
    std::unique_ptr<
        std::map<Vessel const* const,
                 Displacement<Barycentric>>> displacements_from_centre_of_mass;
    std::unique_ptr<
        std::map<Vessel const* const,
                 Velocity<Barycentric>>> velocities_from_centre_of_mass;
    std::unique_ptr<Displacement<World>> displacement_correction;
    std::unique_ptr<Velocity<World>> velocity_correction;
  };

  // Computes the world degrees of freedom of the centre of mass of
  // |next| using the contents of |next->parts|.
  void ComputeNextCentreOfMassWorldDegreesOfFreedom(FullState* next);

  // Computes |next->displacements_from_centre_of_mass| and
  // |next->velocities_from_centre_of_mass|.
  void ComputeNextVesselOffsets(
      PlanetariumRotation const& planetarium_rotation,
      FullState* next);

  // Creates |next->centre_of_mass_trajectory| and appends to it the barycentre
  // of the degrees of freedom of the vessels in |next->vessels|.  There is no
  // intrinsic acceleration.
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
      std::vector<PartCorrespondence>* const common_parts,
      FullState* next);

  // Given the vector of common parts produced by |IntrinsicAcceleration|,
  // constructs |next->centre_of_mass_trajectory| and appends degrees of
  // freedom at |current_time| that conserve the degrees of freedom of the
  // centre of mass of the parts in |common_parts|.  No transfer of ownership.
  void Shift(PlanetariumRotation const& planetarium_rotation,
             Instant const& current_time,
             std::vector<PartCorrespondence> const* const common_parts,
             FullState* next);

  std::unique_ptr<FullState> current_;
  // The following member is only accessed by |AddVesselToNext| and at the
  // beginning of |Prepare|.
  std::unique_ptr<PreliminaryState> next_;

  MasslessBody const body_;
};

}  // namespace ksp_plugin
}  // namespace principia
