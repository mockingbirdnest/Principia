#pragma once

#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/monostable.hpp"
#include "ksp_plugin/vessel.hpp"
#include "ksp_plugin/rendering_frame.hpp"
#include "physics/body.hpp"
#include "physics/n_body_system.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {

using geometry::Displacement;
using geometry::Instant;
using geometry::Point;
using geometry::Rotation;
using integrators::SPRKIntegrator;
using physics::Body;
using physics::Trajectory;
using physics::NBodySystem;
using quantities::Angle;
using si::Second;

// The GUID of a vessel, obtained by |v.id.ToString()| in C#. We use this as a
// key in an |std::map|.
using GUID = std::string;
// The index of a body in |FlightGlobals.Bodies|, obtained by
// |b.flightGlobalsIndex| in C#. We use this as a key in an |std::map|.
using Index = int;

// Represents the line segment {(1-s) |begin| + s |end| | s ∈ [0, 1]}.
// It is immediate that ∀ s ∈ [0, 1], (1-s) |begin| + s |end| is a convex
// combination of |begin| and |end|, so that this is well-defined for |begin|
// and |end| in an affine space.
template<typename Frame>
struct LineSegment {
  LineSegment(Position<Frame> const& begin, Position<Frame> const& end)
      : begin(begin),
        end(end) {}
  Position<Frame> const begin;
  Position<Frame> const end;
};

// We render trajectories as polygons.
template<typename Frame>
using RenderedTrajectory = std::vector<LineSegment<Frame>>;

class Plugin {
 public:
  Plugin() = delete;
  Plugin(Plugin const&) = delete;
  Plugin(Plugin&&) = delete;
  virtual ~Plugin() = default;

  // Constructs a |Plugin|. The current time of that instance is |initial_time|.
  // The angle between the axes of |World| and |Barycentric| at |initial_time|
  // is set to |planetarium_rotation|. Inserts a celestial body with an
  // arbitrary position, index |sun_index| and gravitational parameter
  // |sun_gravitational_parameter|.
  // Starts initialization.
  // The arguments correspond to KSP's
  // |Planetarium.GetUniversalTime()|,
  // |Planetarium.fetch.Sun.flightGlobalsIndex|,
  // |Planetarium.fetch.Sun.gravParameter|,
  // |Planetarium.InverseRotAngle|.
  Plugin(Instant const& initial_time,
         Index const sun_index,
         GravitationalParameter const& sun_gravitational_parameter,
         Angle const& planetarium_rotation);

  // Inserts a new celestial body with index |celestial_index| and gravitational
  // parameter |gravitational_parameter|. No body with index |celestial_index|
  // must already have been inserted. The parent of the new body is the body
  // at index |parent_index|, which must already have been inserted. The state
  // of the new body at current time is given by |AliceSun| offsets from the
  // parent. Must only be called during initialization.
  // For a KSP |CelestialBody| |b|, the arguments correspond to:
  // |b.flightGlobalsIndex|,
  // |b.gravParameter|,
  // |b.orbit.referenceBody.flightGlobalsIndex|,
  // |b.orbit.pos|,
  // |b.orbit.vel|.
  virtual void InsertCelestial(
    Index const celestial_index,
    GravitationalParameter const& gravitational_parameter,
    Index const parent_index,
    Displacement<AliceSun> const& from_parent_position,
    Velocity<AliceSun> const& from_parent_velocity);

  // Ends initialization.
  virtual void EndInitialization();

  // Sets the parent of the celestial body with index |celestial_index| to the
  // one with index |parent_index|. Both bodies must already have been
  // inserted. Must be called after initialization.
  // For a KSP |CelestialBody| |b|, the arguments correspond to
  // |b.flightGlobalsIndex|, |b.orbit.referenceBody.flightGlobalsIndex|.
  virtual void UpdateCelestialHierarchy(Index const celestial_index,
                                        Index const parent_index) const;

  // Inserts a new vessel with GUID |vessel_guid| if it does not already exist,
  // and flags the vessel with GUID |vessel_guid| so it is kept when calling
  // |AdvanceTime|. The parent body for the vessel is set to the one with index
  // |parent_index|. It must already have been inserted using
  // |InsertCelestial|.
  // Returns true if a new vessel was inserted. In that case,
  // |SetVesselStateOffset| must be called with the same GUID before the
  // next call to |AdvanceTime|, |VesselDisplacementFromParent| or
  // |VesselParentRelativeVelocity|, so that the initial state of the new
  // vessel is known. Must be called after initialization.
  // For a KSP |Vessel| |v|, the arguments correspond to
  // |v.id|, |v.orbit.referenceBody.flightGlobalsIndex|.
  virtual bool InsertOrKeepVessel(GUID const& vessel_guid,
                                  Index const parent_index);

  // Set the position and velocity of the vessel with GUID |vessel_guid|
  // relative to its parent at current time. |SetVesselStateOffset| must only
  // be called once per vessel. Must be called after initialization.
  // For a KSP |Vessel| |v|, the arguments correspond to
  // |v.id.ToString()|,
  // |v.orbit.pos|,
  // |v.orbit.vel|.
  virtual void SetVesselStateOffset(
      GUID const& vessel_guid,
      Displacement<AliceSun> const& from_parent_position,
      Velocity<AliceSun> const& from_parent_velocity);

  // Simulates the system until instant |t|. All vessels that have not been
  // refreshed by calling |InsertOrKeepVessel| since the last call to
  // |AdvanceTime| will be removed.  Sets |current_time_| to |t|.
  // Must be called after initialization.  |t| must be greater than
  // |current_time_|.  If |AddVesselToNextPhysicsBubble| was called since the
  // last call to |AdvanceTime|, |BubbleDisplacementCorrection| and
  // |BubbleDisplacementVelocity| must have been called too.
  // |planetarium_rotation| is the value of KSP's |Planetarium.InverseRotAngle|
  // at instant |t|, which provides the rotation between the |World| axes and
  // the |Barycentric| axes (we don't use Planetarium.Rotation since it
  // undergoes truncation to single-precision even though it's a double-
  // precision value).  Note that KSP's |Planetarium.InverseRotAngle| is in
  // degrees.
  virtual void AdvanceTime(Instant const& t, Angle const& planetarium_rotation);

  // Returns the position of the vessel with GUID |vessel_guid| relative to its
  // parent at current time. For a KSP |Vessel| |v|, the argument corresponds to
  // |v.id.ToString()|, the return value to |v.orbit.pos|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept. Must
  // be called after initialization.
  virtual Displacement<AliceSun> VesselDisplacementFromParent(
      GUID const& vessel_guid) const;

  // Returns the velocity of the vessel with GUID |vessel_guid| relative to its
  // parent at current time. For a KSP |Vessel| |v|, the argument corresponds to
  // |v.id.ToString()|, the return value to |v.orbit.vel|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept. Must
  // be called after initialization.
  virtual Velocity<AliceSun> VesselParentRelativeVelocity(
      GUID const& vessel_guid) const;

  // Returns the position of the celestial at index |celestial_index| relative
  // to its parent at current time. For a KSP |CelestialBody| |b|, the argument
  // corresponds to |b.flightGlobalsIndex|, the return value to |b.orbit.pos|.
  // A celestial with index |celestial_index| must have been inserted, and it
  // must not be the sun. Must be called after initialization.
  virtual Displacement<AliceSun> CelestialDisplacementFromParent(
      Index const celestial_index) const;

  // Returns the velocity of the celestial at index |celestial_index| relative
  // to its parent at current time. For a KSP |CelestialBody| |b|, the argument
  // corresponds to |b.flightGlobalsIndex|, the return value to |b.orbit.vel|.
  // A celestial with index |celestial_index| must have been inserted, and it
  // must not be the sun. Must be called after initialization.
  virtual Velocity<AliceSun> CelestialParentRelativeVelocity(
      Index const celestial_index) const;

  // Returns a polygon in |World| space depicting the trajectory of the vessel
  // with the given |GUID| in |frame|.  |sun_world_position| is the current
  // position of the sun in |World| space as returned by
  // |Planetarium.fetch.Sun.position|.  It is used to define the relation
  // between |WorldSun| and |World|.
  virtual RenderedTrajectory<World> RenderedVesselTrajectory(
      GUID const& vessel_guid,
      RenderingFrame const& frame,
      Position<World> const& sun_world_position) const;

  virtual std::unique_ptr<BodyCentredNonRotatingFrame>
  NewBodyCentredNonRotatingFrame(Index const reference_body_index) const;

  virtual std::unique_ptr<BarycentricRotatingFrame>
  NewBarycentricRotatingFrame(Index const primary_index,
                              Index const secondary_index) const;

  virtual Position<World> VesselWorldPosition(
      GUID const& vessel_guid,
      Position<World> const& parent_world_position) const;

  virtual Velocity<World> VesselWorldVelocity(
      GUID const& vessel_guid,
      Velocity<World> const& parent_world_velocity,
      Time const& parent_rotation_period) const;

  // Creates |next_physics_bubble_| if it is null.  Adds the vessel with GUID
  // |vessel_guid| to |next_physics_bubble_->vessels| with a list of pointers to
  // the |Part|s in |parts|.  Merges |parts| into |next_physics_bubble_->parts|.
  // Adds the vessel to |dirty_vessels_|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept.  The
  // vessel with GUID |vessel_guid| must not already be in
  // |next_physics_bubble_->vessels|.  |parts| must not contain a |PartID|
  // already in |next_physics_bubble_->parts|.
  void AddVesselToNextPhysicsBubble(
      GUID const& vessel_guid,
      std::vector<std::pair<PartID, std::unique_ptr<Part<World>>>> parts);
  // Computes and returns |current_physics_bubble_->displacement_correction|.
  // This is the shift to be applied to the physics bubble in worldspace in
  // order for it to be in the correct position.
  Displacement<World> BubbleDisplacementCorrection(
      Position<World> const& sun_world_position) const;
  // Computes and returns |current_physics_bubble_->velocity_correction|.
  // This is the shift to be applied to the physics bubble in worldspace in
  // order for it to be in the correct position.
  Velocity<World> BubbleVelocityCorrection(
      Index const reference_body_index) const;

 private:
  using GUIDToOwnedVessel = std::map<GUID, std::unique_ptr<Vessel>>;
  using GUIDToUnownedVessel = std::map<GUID, Vessel* const>;

  // The common last time of the histories of synchronized vessels and
  // celestials.
  Instant const& HistoryTime() const;

  // The rotation between the |World| basis at |current_time_| and the
  // |Barycentric| axes. Since |WorldSun| is not a rotating reference frame,
  // this change of basis is all that's required to convert relative velocities
  // or displacements between simultaneous events.
  Rotation<Barycentric, WorldSun> PlanetariumRotation() const;

  // Utilities for |AdvanceTime|.

  // Remove vessels not in |kept_|, and clears |kept_|.
  void CleanUpVessels();
  // Given an iterator to an element of |vessels_|, check that the corresponding
  // |Vessel| |is_initialized()|, and that it is not in |new_vessels_| if, and
  // only if, it |is_synchronized()|.
  // Also checks that its |prolongation().last().time()| is at least
  // |HistoryTime()|, and that if it |is_synchronized()|, its
  // |history().last().time()| is exactly |HistoryTime()|.
  void CheckVesselInvariants(GUIDToOwnedVessel::const_iterator const it) const;
  // Evolves the histories of the |celestials_| and of the synchronized vessels
  // up to at most |t|. |t| must be large enough that at least one step of
  // size |Δt_| can fit between |current_time_| and |t|.
  void EvolveHistories(Instant const& t);
  // Synchronizes the |new_vessels_|, clears |new_vessels_|.  Prolongs the
  // histories of the vessels in the physics bubble by evolving the trajectory
  // of the |current_physics_bubble_| if there is one, prolongs the histories of
  // the remaining |dirty_vessels_| using their prolongations, clears
  // |dirty_vessels_|.
  void SynchronizeNewVesselsAndCleanDirtyVessels();
  // Called from |SynchronizeNewVesselsAndCleanDirtyVessels()|, prolongs the
  // histories of the vessels in the physics bubble (the integration must
  // already have been done).  Any new vessels in the physics bubble are
  // synchronized and removed from |new_vessels_|.
  void SynchronizeBubbleHistories();
  // Resets the prolongations of all vessels and celestials to |HistoryTime()|.
  // All vessels must satisfy |is_synchronized()|.
  void ResetProlongations();
  // Evolves the prolongations of all celestials and vessels up to exactly
  // instant |t|.  Also evolves the trajectory/ of the |current_physics_bubble_|
  // if there is one.
  void EvolveProlongationsAndBubble(Instant const& t);
  // Returns true if, and only if, |vessel| is in
  // |current_physics_bubble_->vessels|.  |current_physics_bubble_| may be null,
  // in that case, returns false.
  bool IsInPhysicsBubble(Vessel* const vessel) const;
  // |current_physics_bubble_->vessels.size()|, or 0 if
  // |current_physics_bubble_| is null.
  std::size_t NumberOfVesselsInPhysicsBubble() const;
  // returns |current_physics_bubble_ != nullptr|.
  bool HavePhysicsBubble() const;

  // If |next_physics_bubble_| is not null, computes the world centre of mass,
  // trajectory (including intrinsic acceleration) of |*next_physics_bubble_|.
  // Moves |next_physics_bubble_| into |current_physics_bubble_|.
  void PreparePhysicsBubble(Instant const& next_time);

  // Utilities for |PreparePhysicsBubble|.

  // Computes the world degrees of freedom of the centre of mass of
  // |next_physics_bubble_| using the contents of |next_physics_bubble_->parts|.
  // |next_physics_bubble_| must not be null.
  void ComputeNextPhysicsBubbleCentreOfMassWorldDegreesOfFreedom();
  // Computes |next_physics_bubble_->displacements_from_centre_of_mass| and
  // |next_physics_bubble_->velocities_from_centre_of_mass|.
  // |next_physics_bubble_| must not be null.
  void ComputeNextPhysicsBubbleVesselOffsets();
  // Creates |next_physics_bubble_->centre_of_mass_trajectory| and appends to it
  // the barycentre of the degrees of freedom of the vessels in
  // |next_physics_bubble_->vessels|.  There is no intrinsic acceleration.
  // |next_physics_bubble_| must not be null.
  void RestartNextPhysicsBubble();
  // Returns the intrinsic acceleration measured on the parts that are common to
  // the current and next physics bubbles.  Stores a pair of pointers to parts
  // (current, next) in |common_parts| for all parts common to the current and
  // next physics bubble.
  // |common_parts| must not be null.  |*common_parts| must be empty.  No
  // transfer of ownership.
  Vector<Acceleration, World> IntrinsicAcceleration(
      Instant const& next_time,
      std::vector<std::pair<Part<World>*, Part<World>*>>* const common_parts);
  // Given the vector of common parts produced by |IntrinsicAcceleration|,
  // constructs |*next_physics_bubble_->centre_of_mass_trajectory| and appends
  // degres of freedom at |current_time_| that conserve the degrees of freedom
  // of the centre of mass of the parts in |common_parts|.
  // |common_parts| must not be null.  |next_physics_bubble_| must not be null.
  // No transfer of ownership.
  void ShiftBubble(
      std::vector<std::pair<Part<World>*,
                            Part<World>*>> const* const common_parts);

  // TODO(egg): Constant time step for now.
  Time const Δt_ = 10 * Second;

  GUIDToOwnedVessel vessels_;
  std::map<Index, std::unique_ptr<Celestial>> celestials_;

  // The vessels which have been inserted after |HistoryTime()|.  These are the
  // vessels which do not satisfy |is_synchronized()|, i.e., they do not have a
  // history.  The pointers are not owning and not null.
  std::set<Vessel* const> new_vessels_;
  // The vessels that have been added to the physics bubble after
  // |HistoryTime()|.  For these vessels, the prolongation contains information
  // that may not be discarded, and the history will be advanced using the
  // prolongation.  The pointers are not owning and not null.
  std::set<Vessel* const> dirty_vessels_;

  // The vessels that will be kept during the next call to |AdvanceTime|.
  std::set<Vessel const* const> kept_;

  struct PhysicsBubble {
    std::map<Vessel* const, std::vector<Part<World>* const>> vessels;
    std::map<PartID, std::unique_ptr<Part<World>> const> parts;
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
    std::unique_ptr<Displacement<World>> displacement_correction;
    std::unique_ptr<Velocity<World>> velocity_correction;
  };

  std::unique_ptr<PhysicsBubble> current_physics_bubble_;
  std::unique_ptr<PhysicsBubble> next_physics_bubble_;

  MasslessBody bubble_body_;

  std::unique_ptr<NBodySystem<Barycentric>> n_body_system_;
  // The symplectic integrator computing the synchronized histories.
  SPRKIntegrator<Length, Speed> history_integrator_;
  // The integrator computing the prolongations.
  SPRKIntegrator<Length, Speed> prolongation_integrator_;

  // Whether initialization is ongoing.
  Monostable initializing;

  Angle planetarium_rotation_;
  // The current in-game universal time.
  Instant current_time_;
  Celestial* sun_;  // Not owning, not null.

  friend class TestablePlugin;
};

}  // namespace ksp_plugin
}  // namespace principia
