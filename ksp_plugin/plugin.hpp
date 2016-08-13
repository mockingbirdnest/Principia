
#pragma once

#include <limits>
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
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/physics_bubble.hpp"
#include "ksp_plugin/vessel.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/frame_field.hpp"
#include "physics/hierarchical_system.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using geometry::Point;
using geometry::Rotation;
using integrators::FixedStepSizeIntegrator;
using integrators::AdaptiveStepSizeIntegrator;
using physics::Body;
using physics::DiscreteTrajectory;
using physics::DynamicFrame;
using physics::Ephemeris;
using physics::FrameField;
using physics::Frenet;
using physics::HierarchicalSystem;
using physics::RelativeDegreesOfFreedom;
using quantities::Angle;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;

namespace ksp_plugin {

// The GUID of a vessel, obtained by |v.id.ToString()| in C#. We use this as a
// key in an |std::map|.
using GUID = std::string;
// The index of a body in |FlightGlobals.Bodies|, obtained by
// |b.flightGlobalsIndex| in C#. We use this as a key in an |std::map|.
using Index = int;

class Plugin {
 public:
  Plugin() = delete;
  Plugin(Plugin const&) = delete;
  Plugin(Plugin&&) = delete;
  Plugin& operator=(Plugin const&) = delete;
  Plugin& operator=(Plugin&&) = delete;
  virtual ~Plugin() = default;

  // Constructs a |Plugin|. The current time of that instance is
  // |solar_system_epoch|.  The angle between the axes of |World| and
  // |Barycentric| at |solar_system_epoch| is set to |planetarium_rotation|.
  Plugin(Instant const& game_epoch,
         Instant const& solar_system_epoch,
         Angle const& planetarium_rotation);

  // Inserts a celestial body with index |celestial_index| body |body|,
  // giving it the initial state |initial_state|.
  // If |parent_index| is null, inserts the sun, otherwise the parent of the new
  // body is the body with index |*parent_index|, which must already have been
  // inserted.  Hierarchical initialization must not be ongoing.
  virtual void InsertCelestialAbsoluteCartesian(
      Index const celestial_index,
      std::experimental::optional<Index> const& parent_index,
      DegreesOfFreedom<Barycentric> const& initial_state,
      not_null<std::unique_ptr<MassiveBody const>> body);

  // Hierarchical initialization must be ongoing.
  virtual void InsertCelestialJacobiKeplerian(
      Index const celestial_index,
      std::experimental::optional<Index> const& parent_index,
      std::experimental::optional<
          physics::KeplerianElements<Barycentric>> const& keplerian_elements,
      not_null<std::unique_ptr<MassiveBody>> body);

  // Ends initialization.  The sun must have been inserted.
  virtual void EndInitialization();

  // Returns true iff this is the unstable KSP stock system.  Must be called
  // after initialization.
  virtual bool IsKspStockSystem() const;

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

  // Forgets the histories of the |celestials_| and of the vessels before |t|.
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

  // Updates the prediction for the vessel with guid |vessel_guid|.
  void UpdatePrediction(GUID const& vessel_guid) const;

  virtual void CreateFlightPlan(GUID const& vessel_guid,
                                Instant const& final_time,
                                Mass const& initial_mass) const;

  // Returns a polygon in |World| space depicting the trajectory of the vessel
  // with the given |GUID| in the frame defined by the current
  // |plotting_frame_|.
  // |sun_world_position| is the current position of the sun in |World| space as
  // returned by |Planetarium.fetch.Sun.position|.  It is used to define the
  // relation between |WorldSun| and |World|.  No transfer of ownership.
  virtual not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderedVesselTrajectory(GUID const& vessel_guid,
                           Position<World> const& sun_world_position) const;

  // Returns a polygon in |World| space depicting the trajectory of
  // |predicted_vessel_| from |CurrentTime()| to
  // |CurrentTime() + prediction_length_| in the current |plotting_frame_|.
  // |sun_world_position| is the current position of the sun in |World| space as
  // returned by |Planetarium.fetch.Sun.position|.  It is used to define the
  // relation between |WorldSun| and |World|.
  // No transfer of ownership.
  // |predicted_vessel_| must have been set, and |AdvanceTime()| must have been
  // called after |predicted_vessel_| was set.
  virtual not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderedPrediction(GUID const& vessel_guid,
                     Position<World> const& sun_world_position) const;

  // A utility for |RenderedPrediction| and |RenderedVesselTrajectory|,
  // returns a |Positions| object corresponding to the trajectory defined by
  // |begin| and |end|, as seen in the current |plotting_frame_|.
  // TODO(phl): Use this directly in the interface and remove the other
  // |Rendered...|.
  virtual not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderedTrajectoryFromIterators(
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Position<World> const& sun_world_position) const;

  virtual void ComputeAndRenderApsides(
      Index const celestial_index,
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Position<World> const& sun_world_position,
      std::unique_ptr<DiscreteTrajectory<World>>& apoapsides,
      std::unique_ptr<DiscreteTrajectory<World>>& periapsides) const;

  virtual void SetPredictionLength(Time const& t);

  virtual void SetPredictionAdaptiveStepParameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          prediction_adaptive_step_parameters);

  virtual bool HasVessel(GUID const& vessel_guid) const;
  virtual not_null<Vessel*> GetVessel(GUID const& vessel_guid) const;

  virtual not_null<std::unique_ptr<NavigationFrame>>
  NewBodyCentredNonRotatingNavigationFrame(
      Index const reference_body_index) const;

  virtual not_null<std::unique_ptr<NavigationFrame>>
  NewBarycentricRotatingNavigationFrame(Index const primary_index,
                                        Index const secondary_index) const;

  virtual void SetPlottingFrame(
      not_null<std::unique_ptr<NavigationFrame>> plotting_frame);
  virtual not_null<NavigationFrame const*> GetPlottingFrame() const;

  virtual Position<World> PlotBarycentricPosition(
      Instant const& t,
      Position<Barycentric> const& position,
      Position<World> const& sun_world_position) const;

  // Creates |next_physics_bubble_| if it is null.  Adds the vessel with GUID
  // |vessel_guid| to |next_physics_bubble_->vessels| with a list of pointers to
  // the |Part|s in |parts|.  Merges |parts| into |next_physics_bubble_->parts|.
  // Marks the vessel as dirty.
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

  // The navball field at |current_time| for the current |plotting_frame_|.
  virtual FrameField<World> Navball(
      Position<World> const& sun_world_position) const;

  // The unit tangent, normal, or binormal vector to the trajectory of the
  // vessel with the given GUID in the current |plotting_frame_|.
  virtual Vector<double, World> VesselTangent(GUID const& vessel_guid) const;
  virtual Vector<double, World> VesselNormal(GUID const& vessel_guid) const;
  virtual Vector<double, World> VesselBinormal(GUID const& vessel_guid) const;

  virtual Velocity<World> VesselVelocity(GUID const& vessel_guid) const;

  // Returns
  // |sun_looking_glass.Inverse().Forget() * PlanetariumRotation().Forget()|.
  virtual OrthogonalMap<Barycentric, WorldSun> BarycentricToWorldSun() const;

  virtual Instant GameEpoch() const;

  virtual Instant CurrentTime() const;

  // Must be called after initialization.
  virtual void WriteToMessage(
      not_null<serialization::Plugin*> const message) const;
  static not_null<std::unique_ptr<Plugin>> ReadFromMessage(
      serialization::Plugin const& message);

 private:
  using GUIDToOwnedVessel = std::map<GUID, not_null<std::unique_ptr<Vessel>>>;
  using GUIDToUnownedVessel = std::map<GUID, not_null<Vessel*> const>;
  using IndexToOwnedCelestial =
      std::map<Index, not_null<std::unique_ptr<Celestial>>>;
  using NewtonianMotionEquation =
      Ephemeris<Barycentric>::NewtonianMotionEquation;
  using IndexToMassiveBody =
      std::map<Index, not_null<std::unique_ptr<MassiveBody const>>>;
  using IndexToDegreesOfFreedom =
      std::map<Index, DegreesOfFreedom<Barycentric>>;
  using Trajectories = std::vector<not_null<DiscreteTrajectory<Barycentric>*>>;

  // This constructor should only be used during deserialization.  All vessels
  // are added to |kept_vessels_|  The resulting plugin is not |initializing_|.
  Plugin(GUIDToOwnedVessel vessels,
         IndexToOwnedCelestial celestials,
         not_null<std::unique_ptr<PhysicsBubble>> bubble,
         std::unique_ptr<Ephemeris<Barycentric>> ephemeris,
         Ephemeris<Barycentric>::FixedStepParameters const& history_parameters,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prolongation_parameters,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prediction_parameters,
         Angle const& planetarium_rotation,
         Instant const& game_epoch,
         Instant const& current_time,
         Index sun_index);

  // We virtualize this function for testing purposes.
  // Requires |absolute_initialization_| and consumes it.
  virtual void InitializeEphemerisAndSetCelestialTrajectories();

  not_null<std::unique_ptr<Vessel>> const& find_vessel_by_guid_or_die(
      GUID const& vessel_guid) const;

  // The rotation between the |AliceWorld| basis at |current_time_| and the
  // |Barycentric| axes. Since |AliceSun| is not a rotating reference frame,
  // this change of basis is all that's required to convert relative velocities
  // or displacements between simultaneous events.
  Rotation<Barycentric, AliceSun> PlanetariumRotation() const;

  // Utilities for |AdvanceTime|.

  // Remove vessels not in |kept_vessels_|, and clears |kept_vessels_|.
  void FreeVessels();
  // Evolves the trajectory of the |current_physics_bubble_|.
  void EvolveBubble(Instant const& t);

  Vector<double, World> FromVesselFrenetFrame(
      Vessel const& vessel,
      Vector<double, Frenet<Navigation>> const& vector) const;

  // Fill |celestials| using the |index| and |parent_index| fields found in
  // |celestial_messages| (which may be pre- or post-Bourbaki).
  template<typename T>
  static void ReadCelestialsFromMessages(
    Ephemeris<Barycentric> const& ephemeris,
    google::protobuf::RepeatedPtrField<T> const& celestial_messages,
    not_null<IndexToOwnedCelestial*> const celestials);

  // Computes a fingerprint for the parameters passed to
  // |InsertCelestialJacobiKeplerian|.
  static std::uint64_t FingerprintCelestialJacobiKeplerian(
      Index const celestial_index,
      std::experimental::optional<Index> const& parent_index,
      std::experimental::optional<
          physics::KeplerianElements<Barycentric>> const& keplerian_elements,
      MassiveBody const& body);

  GUIDToOwnedVessel vessels_;
  IndexToOwnedCelestial celestials_;

  // The vessels that will be kept during the next call to |AdvanceTime|.
  std::set<not_null<Vessel const*>> kept_vessels_;

  not_null<std::unique_ptr<PhysicsBubble>> const bubble_;

  struct AbsoluteInitializationObjects {
    IndexToMassiveBody bodies;
    IndexToDegreesOfFreedom initial_state;
  };
  std::experimental::optional<AbsoluteInitializationObjects>
      absolute_initialization_;

  struct HierarchicalInitializationObjects {
    HierarchicalInitializationObjects(
        not_null<std::unique_ptr<MassiveBody const>> sun)
        : system(std::move(sun)) {}
    HierarchicalSystem<Barycentric> system;
    std::map<Index, MassiveBody const*> indices_to_bodies;
    std::map<Index, std::experimental::optional<Index>> parents;
  };
  std::experimental::optional<HierarchicalInitializationObjects>
      hierarchical_initialization_;

  // Null if and only if |initializing_|.
  // TODO(egg): optional.
  std::unique_ptr<Ephemeris<Barycentric>> ephemeris_;

  // The parameters for computing the various trajectories.
  Ephemeris<Barycentric>::FixedStepParameters history_parameters_;
  Ephemeris<Barycentric>::AdaptiveStepParameters prolongation_parameters_;
  Ephemeris<Barycentric>::AdaptiveStepParameters prediction_parameters_;
  Time prediction_length_ = 1 * Hour;

  // Whether initialization is ongoing.
  base::Monostable initializing_;

  Angle planetarium_rotation_;
  // The game epoch in real time.
  Instant const game_epoch_;
  // The current in-game universal time.
  Instant current_time_;

  Celestial* sun_ = nullptr;  // Not owning, not null after InsertSun is called.

  // Not null after initialization. |EndInitialization| sets it to the
  // heliocentric frame.
  std::unique_ptr<NavigationFrame> plotting_frame_;

  // Used for detecting and patching the stock system.
  std::set<std::uint64_t> celestial_jacobi_keplerian_fingerprints_;
  bool is_ksp_stock_system_ = false;

  friend class TestablePlugin;
};

}  // namespace ksp_plugin
}  // namespace principia
