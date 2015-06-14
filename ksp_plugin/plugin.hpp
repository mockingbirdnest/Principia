#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "base/monostable.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/physics_bubble.hpp"
#include "ksp_plugin/vessel.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/body.hpp"
#include "physics/ephemeris.hpp"
#include "physics/trajectory.hpp"
#include "physics/transforms.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {

using geometry::Displacement;
using geometry::Instant;
using geometry::Point;
using geometry::Rotation;
using integrators::FixedStepSizeIntegrator;
using integrators::AdaptiveStepSizeIntegrator;
using physics::Body;
using physics::FrameField;
using physics::Ephemeris;
using physics::Trajectory;
using physics::Transforms;
using quantities::Angle;
using si::Hour;
using si::Metre;
using si::Milli;
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

using RenderingTransforms =
    Transforms<MobileInterface, Barycentric, Rendering, Barycentric>;

class Plugin {
 public:
  Plugin() = delete;
  Plugin(Plugin const&) = delete;
  Plugin(Plugin&&) = delete;
  Plugin& operator=(Plugin const&) = delete;
  Plugin& operator=(Plugin&&) = delete;
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
  // |{b.orbit.pos, b.orbit.vel}|.
  virtual void InsertCelestial(
    Index const celestial_index,
    GravitationalParameter const& gravitational_parameter,
    Index const parent_index,
    RelativeDegreesOfFreedom<AliceSun> const& from_parent);

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
  // |{v.orbit.pos, v.orbit.vel}|.
  virtual void SetVesselStateOffset(
      GUID const& vessel_guid,
      RelativeDegreesOfFreedom<AliceSun> const& from_parent);

  // Simulates the system until instant |t|. All vessels that have not been
  // refreshed by calling |InsertOrKeepVessel| since the last call to
  // |AdvanceTime| will be removed.  Sets |current_time_| to |t|.
  // Must be called after initialization.  |t| must be greater than
  // |current_time_|.  If |PhysicsBubble::AddVesselToNext| was called since the
  // last call to |AdvanceTime|, |PhysicsBubble::DisplacementCorrection| and
  // |PhysicsBubble::DisplacementVelocity| must have been called too.
  // |planetarium_rotation| is the value of KSP's |Planetarium.InverseRotAngle|
  // at instant |t|, which provides the rotation between the |World| axes and
  // the |Barycentric| axes (we don't use Planetarium.Rotation since it
  // undergoes truncation to single-precision even though it's a double-
  // precision value).  Note that KSP's |Planetarium.InverseRotAngle| is in
  // degrees.
  virtual void AdvanceTime(Instant const& t, Angle const& planetarium_rotation);

  // Forgets the histories of the |celestials_| and of the synchronized vessels
  // before |t|.
  virtual void ForgetAllHistoriesBefore(Instant const& t) const;

  // Returns the displacement and velocity of the vessel with GUID |vessel_guid|
  // relative to its parent at current time. For a KSP |Vessel| |v|, the
  // argument corresponds to  |v.id.ToString()|, the return value to
  // |{v.orbit.pos, v.orbit.vel}|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept. Must
  // be called after initialization.
  virtual RelativeDegreesOfFreedom<AliceSun> VesselFromParent(
      GUID const& vessel_guid) const;

  // Returns the displacement and velocity of the celestial at index
  // |celestial_index| relative to its parent at current time. For a KSP
  // |CelestialBody| |b|, the argument corresponds to |b.flightGlobalsIndex|,
  // the return value to |{b.orbit.pos, b.orbit.vel}|.
  // A celestial with index |celestial_index| must have been inserted, and it
  // must not be the sun. Must be called after initialization.
  virtual RelativeDegreesOfFreedom<AliceSun> CelestialFromParent(
      Index const celestial_index) const;

  // Returns a polygon in |World| space depicting the trajectory of the vessel
  // with the given |GUID| in the frame defined by |transforms|.
  // |sun_world_position| is the current position of the sun in |World| space as
  // returned by |Planetarium.fetch.Sun.position|.  It is used to define the
  // relation between |WorldSun| and |World|.  No transfer of ownership.
  virtual RenderedTrajectory<World> RenderedVesselTrajectory(
      GUID const& vessel_guid,
      not_null<RenderingTransforms*> const transforms,
      Position<World> const& sun_world_position) const;

  // Returns a polygon in |World| space depicting the trajectory of
  // |predicted_vessel_| from |current_time()| to
  // |current_time() + prediction_length_| in the frame defined by |transforms|.
  // |sun_world_position| is the current  position of the sun in |World| space
  // as returned by |Planetarium.fetch.Sun.position|.  It is used to define the
  // relation between |WorldSun| and |World|.  No transfer of ownership.
  // |predicted_vessel_| must have been set, and |AdvanceTime()| must have been
  // called after |predicted_vessel_| was set.  Not const because of the stupid
  // global variable |transforms_are_operating_on_predictions_|.
  virtual RenderedTrajectory<World> RenderedPrediction(
      not_null<RenderingTransforms*> const transforms,
      Position<World> const& sun_world_position);

  virtual void set_predicted_vessel(GUID const& vessel_guid);
  // Calls |DeletePredictions()| and nulls |predicted_vessel_|.
  virtual void clear_predicted_vessel();

  virtual void set_prediction_length(Time const& t);

  // The step used when computing the prediction.
  virtual void set_prediction_step(Time const& t);

  virtual bool has_vessel(GUID const& vessel_guid) const;

  virtual not_null<std::unique_ptr<RenderingTransforms>>
  NewBodyCentredNonRotatingTransforms(Index const reference_body_index) const;

  virtual not_null<std::unique_ptr<RenderingTransforms>>
  NewBarycentricRotatingTransforms(Index const primary_index,
                                   Index const secondary_index) const;

  // Creates |next_physics_bubble_| if it is null.  Adds the vessel with GUID
  // |vessel_guid| to |next_physics_bubble_->vessels| with a list of pointers to
  // the |Part|s in |parts|.  Merges |parts| into |next_physics_bubble_->parts|.
  // Adds the vessel to |dirty_vessels_|.
  // A vessel with GUID |vessel_guid| must have been inserted and kept.  The
  // vessel with GUID |vessel_guid| must not already be in
  // |next_physics_bubble_->vessels|.  |parts| must not contain a |PartId|
  // already in |next_physics_bubble_->parts|.
  virtual void AddVesselToNextPhysicsBubble(GUID const& vessel_guid,
                                            std::vector<IdAndOwnedPart> parts);

  // Returns |bubble_.empty()|.
  virtual bool PhysicsBubbleIsEmpty() const;

  // Computes and returns |current_physics_bubble_->displacement_correction|.
  // This is the |World| shift to be applied to the physics bubble in order for
  // it to be in the correct position.
  virtual Displacement<World> BubbleDisplacementCorrection(
      Position<World> const& sun_world_position) const;
  // Computes and returns |current_physics_bubble_->velocity_correction|.
  // This is the |World| shift to be applied to the physics bubble in order for
  // it to have the correct velocity.
  virtual Velocity<World> BubbleVelocityCorrection(
      Index const reference_body_index) const;

  // The navball field at |current_time| for the given |transforms|.
  virtual FrameField<World> Navball(
      not_null<RenderingTransforms*> const transforms,
      Position<World> const& sun_world_position) const;

  // The unit tangent vector to the trajectory of the vessel with the given GUID
  // in the frame given by |transforms|.
  virtual Vector<double, World> VesselTangent(
      GUID const& vessel_guid,
      not_null<RenderingTransforms*> const transforms) const;

  virtual Instant current_time() const;

  // Must be called after initialization.
  virtual void WriteToMessage(
      not_null<serialization::Plugin*> const message) const;
  // NOTE(egg): This should return a |not_null|, but we can't do that until
  // |not_null<std::unique_ptr<T>>| is convertible to |std::unique_ptr<T>|, and
  // that requires a VS 2015 feature (rvalue references for |*this|).
  static std::unique_ptr<Plugin> ReadFromMessage(
      serialization::Plugin const& message);

 private:
  using GUIDToOwnedVessel = std::map<GUID, not_null<std::unique_ptr<Vessel>>>;
  using GUIDToUnownedVessel = std::map<GUID, not_null<Vessel*> const>;
  using IndexToOwnedCelestial =
      std::map<Index, not_null<std::unique_ptr<Celestial>>>;
  using NewtonianMotionEquation =
      Ephemeris<Barycentric>::NewtonianMotionEquation;
  using CelestialToMassiveBody =
      std::map<Celestial const*, std::unique_ptr<MassiveBody const>>;
  using CelestialToDegreesOfFreedom =
      std::map<Celestial const*, DegreesOfFreedom<Barycentric>>;

  // This constructor should only be used during deserialization.
  // |unsynchronized_vessels_| is initialized consistently.  All vessels are
  // added to |kept_vessels_|  The resulting plugin is not |initializing_|.
  Plugin(GUIDToOwnedVessel vessels,
         IndexToOwnedCelestial celestials,
         std::set<not_null<Vessel*> const> dirty_vessels,
         not_null<std::unique_ptr<PhysicsBubble>> bubble,
         Angle planetarium_rotation,
         Instant current_time,
         Index sun_index);

  not_null<std::unique_ptr<Vessel>> const& find_vessel_by_guid_or_die(
      GUID const& vessel_guid) const;

  // Returns |!dirty_vessels_.empty()|.
  bool has_dirty_vessels() const;
  // Returns |!unsynchronized_vessels_.empty()|.
  bool has_unsynchronized_vessels() const;
  // Returns |dirty_vessels_.count(vessel) > 0|.
  bool is_dirty(not_null<Vessel*> const vessel) const;
  // Returns |predicted_vessel_ != nullptr|.
  bool has_predicted_vessel() const;
  // Returns true if there is a |predicted_vessel_| and it has a prediction.
  bool HasPredictions() const;
  // Deletes all the predictions.
  void DeletePredictions();

  // The common last time of the histories of synchronized vessels and
  // celestials.
  Instant const& HistoryTime() const;

  // The rotation between the |AliceWorld| basis at |current_time_| and the
  // |Barycentric| axes. Since |AliceSun| is not a rotating reference frame,
  // this change of basis is all that's required to convert relative velocities
  // or displacements between simultaneous events.
  Rotation<Barycentric, AliceSun> PlanetariumRotation() const;

  // returns
  // |kSunLookingGlass.Inverse().Forget() * PlanetariumRotation().Forget()|.
  OrthogonalMap<Barycentric, WorldSun> BarycentricToWorldSun() const;

  // Utilities for |AdvanceTime|.

  // Remove vessels not in |kept_vessels_|, and clears |kept_vessels_|.
  void CleanUpVessels();
  // Given an iterator to an element of |vessels_|, check that the corresponding
  // |Vessel| |is_initialized()|, and that it is not in
  // |unsynchronized_vessels_| if, and only if, it |is_synchronized()|.
  // Also checks that its |prolongation().last().time()| is at least
  // |HistoryTime()|, and that if it |is_synchronized()|, its
  // |history().last().time()| is exactly |HistoryTime()|.
  void CheckVesselInvariants(GUIDToOwnedVessel::const_iterator const it) const;
  // Evolves the histories of the |celestials_| and of the synchronized vessels
  // up to at most |t|. |t| must be large enough that at least one step of
  // size |Δt_| can fit between |current_time_| and |t|.
  void EvolveHistories(Instant const& t);
  // Synchronizes the |unsynchronized_vessels_|, clears
  // |unsynchronized_vessels_|.  Prolongs the histories of the vessels in the
  // physics bubble by evolving the trajectory of the |current_physics_bubble_|
  // if there is one, prolongs the histories of the remaining |dirty_vessels_|
  // using their prolongations, clears |dirty_vessels_|.
  void SynchronizeNewVesselsAndCleanDirtyVessels();
  // Called from |SynchronizeNewVesselsAndCleanDirtyVessels()|, prolongs the
  // histories of the vessels in the physics bubble (the integration must
  // already have been done).  Any new vessels in the physics bubble are
  // synchronized and removed from |unsynchronized_vessels_|.
  void SynchronizeBubbleHistories();
  // Resets the prolongations of all vessels and celestials to |HistoryTime()|.
  // All vessels must satisfy |is_synchronized()|.
  void ResetProlongations();
  // Evolves the prolongations of all celestials and vessels up to exactly
  // instant |t|.  Also evolves the trajectory of the |current_physics_bubble_|
  // if there is one.
  void EvolveProlongationsAndBubble(Instant const& t);
  // Calls |DeletePredictions()|.  If |has_predicted_vessel()|, computes
  // |system_predictions_| and |prediction_| for the |predicted_vessel_|
  // according to |prediction_length_| and |prediction_step_|.
  void UpdatePredictions();

  // A utility for |RenderedPrediction| and |RenderedVesselTrajectory|,
  // returns a |RenderedTrajectory| as computed by the given |transforms|
  // from the trajectory of |body| starting at |actual_it|.
  RenderedTrajectory<World> RenderTrajectory(
      not_null<Body const*> const body,
      Trajectory<Barycentric>::TransformingIterator<Rendering> const& actual_it,
      not_null<RenderingTransforms*>const transforms,
      Position<World> const& sun_world_position) const;

  Time const Δt_ = 10 * Second;
  Length const prolongation_length_tolerance = 1 * Milli(Metre);
  Speed const prolongation_speed_tolerance = 1 * Milli(Metre) / Second;

  GUIDToOwnedVessel vessels_;
  IndexToOwnedCelestial celestials_;

  // The vessels which have been inserted after |HistoryTime()|.  These are the
  // vessels which do not satisfy |is_synchronized()|, i.e., they do not have a
  // history.  The pointers are not owning.
  std::set<not_null<Vessel*> const> unsynchronized_vessels_;
  // The vessels that have been added to the physics bubble after
  // |HistoryTime()|.  For these vessels, the prolongation contains information
  // that may not be discarded, and the history will be advanced using the
  // prolongation.  The pointers are not owning.
  std::set<not_null<Vessel*> const> dirty_vessels_;

  // The vessels that will be kept during the next call to |AdvanceTime|.
  std::set<not_null<Vessel const*> const> kept_vessels_;

  // Only one prediction for now, using constant timestep.
  Vessel* predicted_vessel_ = nullptr;
  Time prediction_length_ = 1 * Hour;
  Length prediction_length_tolerance = 1 * Metre;
  Speed prediction_speed_tolerance = 1 * Metre / Second;

  not_null<std::unique_ptr<PhysicsBubble>> const bubble_;

  // |bodies_| and |initial_state_| are null if and only if |!initializing_|.
  // TODO(egg): optional.
  std::unique_ptr<CelestialToMassiveBody> bodies_;
  std::unique_ptr<CelestialToDegreesOfFreedom> initial_state_;
  // Null if and only if |initializing_|.
  // TODO(egg): optional.
  std::unique_ptr<Ephemeris<Barycentric>> n_body_system_;
  // The integrator computing the synchronized histories of the vessels.
  FixedStepSizeIntegrator<NewtonianMotionEquation> const& history_integrator_;
  // The integrator computing the prolongations.
  AdaptiveStepSizeIntegrator<
      NewtonianMotionEquation> const& prolongation_integrator_;
  // The integrator computing the predictions.
  AdaptiveStepSizeIntegrator<
      NewtonianMotionEquation> const& prediction_integrator_;

  // Whether initialization is ongoing.
  base::Monostable initializing_;

  Angle planetarium_rotation_;
  // The current in-game universal time.
  Instant current_time_;

  Celestial* sun_;  // Not owning, not null after construction.

  friend class TestablePlugin;
};

}  // namespace ksp_plugin
}  // namespace principia
