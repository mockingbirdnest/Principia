
#pragma once

#include <limits>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "base/monostable.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
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
namespace ksp_plugin {
namespace internal_plugin {

using base::not_null;
using base::Subset;
using geometry::AffineMap;
using geometry::Displacement;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Point;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using integrators::FixedStepSizeIntegrator;
using integrators::AdaptiveStepSizeIntegrator;
using physics::Body;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::DynamicFrame;
using physics::Ephemeris;
using physics::FrameField;
using physics::Frenet;
using physics::HierarchicalSystem;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using physics::RotatingBody;
using quantities::Angle;
using quantities::Force;
using quantities::Length;
using quantities::Mass;
using quantities::Time;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;

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
  virtual ~Plugin();

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
      Index celestial_index,
      std::experimental::optional<Index> const& parent_index,
      DegreesOfFreedom<Barycentric> const& initial_state,
      not_null<std::unique_ptr<MassiveBody const>> body);

  // Hierarchical initialization must be ongoing.
  virtual void InsertCelestialJacobiKeplerian(
      Index celestial_index,
      std::experimental::optional<Index> const& parent_index,
      std::experimental::optional<
          physics::KeplerianElements<Barycentric>> const& keplerian_elements,
      not_null<std::unique_ptr<MassiveBody>> body);

  // Ends initialization.  The sun must have been inserted.
  virtual void EndInitialization();

  // Returns true iff this is the unstable KSP stock system.  Must be called
  // after initialization.
  virtual bool IsKspStockSystem() const;

  // And there shall, in that time, be rumors of things going astray, and there
  // will be a great confusion as to where things really are, and nobody will
  // really know where lieth those little things with the sort of raffia work
  // base, that has an attachment.
  virtual bool HasEncounteredApocalypse(std::string* details) const;

  // Sets the parent of the celestial body with index |celestial_index| to the
  // one with index |parent_index|. Both bodies must already have been
  // inserted. Must be called after initialization.
  // For a KSP |CelestialBody| |b|, the arguments correspond to
  // |b.flightGlobalsIndex|, |b.orbit.referenceBody.flightGlobalsIndex|.
  virtual void UpdateCelestialHierarchy(Index celestial_index,
                                        Index parent_index) const;

  // Sets the celestial whose axis of rotation will coincide with the |Alice|
  // z axis.
  virtual void SetMainBody(Index index);
  virtual Rotation<BodyWorld, World> CelestialRotation(Index index) const;
  virtual Rotation<CelestialSphere, World> CelestialSphereRotation() const;

  virtual Angle CelestialInitialRotation(Index celestial_index) const;
  virtual Time CelestialRotationPeriod(Index celestial_index) const;

  // Inserts a new vessel with GUID |vessel_guid| if it does not already exist,
  // and flags the vessel with GUID |vessel_guid| so it is kept when calling
  // |FreeVesselsAndPartsAndCollectPileUps|. The parent body for the vessel is
  // set to the one with index |parent_index|, which must have been inserted
  // during initialization.
  // Sets |inserted| to true if a new vessel was inserted, to false otherwise.
  // If |InsertOrKeepVessel| is called with |loaded=false|, and returns
  // |inserted=true|, |InsertUnloadedPart| must be called for its parts
  // before the call to |AdvanceTime|, giving the vessel an initial state.
  // For a KSP |Vessel| |v|, the arguments correspond to |v.id|,
  // |v.orbit.referenceBody.flightGlobalsIndex|, |v.loaded|.
  virtual void InsertOrKeepVessel(GUID const& vessel_guid,
                                  std::string const& vessel_name,
                                  Index parent_index,
                                  bool loaded,
                                  bool& inserted);

  // Adds a part with the given |part_id| to the vessel with the given |GUID|,
  // which must be unloaded, putting the part at the given offset from the
  // parent body of the vessel.  The part is given unit mass; this does not
  // matter, since the |PileUp| will be deformed when it first loads anyway.
  virtual void InsertUnloadedPart(
      PartId part_id,
      std::string const& name,
      GUID const& vessel_guid,
      RelativeDegreesOfFreedom<AliceSun> const& from_parent);

  // Inserts a new part with the given ID if it does not already exist, and
  // flags the part with the given ID so it is kept when calling
  // |FreeVesselsAndPartsAndCollectPileUps|.
  // The part is created in the given |vessel|, or if it already existed in
  // another vessel, is moved to that one.  If the part is created, its degrees
  // of freedom are set using those given and the |main_body_index|; otherwise
  // these three parameters are ignored.
  virtual void InsertOrKeepLoadedPart(
      PartId part_id,
      std::string const& name,
      Mass const& mass,
      GUID const& vessel_guid,
      Index main_body_index,
      DegreesOfFreedom<World> const& main_body_degrees_of_freedom,
      DegreesOfFreedom<World> const& part_degrees_of_freedom);

  // Calls |increment_intrinsic_force| on the relevant part, which must be in a
  // loaded vessel.
  virtual void IncrementPartIntrinsicForce(PartId part_id,
                                           Vector<Force, World> const& force);

  // Calls |MakeSingleton| for all parts in loaded vessels, enabling the use of
  // union-find for pile up construction.  This must be called after the calls
  // to |IncrementPartIntrinsicForce|, and before the calls to
  // |ReportCollision|.
  virtual void PrepareToReportCollisions();

  // Notifies |this| that the given vessels are touching, and should gravitate
  // as part of a single rigid body.
  virtual void ReportCollision(PartId part1, PartId part2) const;

  // Destroys the vessels for which |InsertOrKeepVessel| has not been called
  // since the last call to |FreeVesselsAndCollectPileUps|, as well as the parts
  // in loaded vessels for which |InsertOrKeepLoadedPart| has not been called,
  // and updates the list of |pile_ups_| according to the reported collisions.
  virtual void FreeVesselsAndPartsAndCollectPileUps();

  // Calls |SetPartApparentDegreesOfFreedom| on the pile-up containing the
  // relevant part.  This part must be in a loaded vessel.
  virtual void SetPartApparentDegreesOfFreedom(
      PartId part_id,
      DegreesOfFreedom<World> const& degrees_of_freedom);

  // Simulates the system until instant |t|.  Sets |current_time_| to |t|.
  // Must be called after initialization.
  // Clears the intrinsic force on all loaded parts.
  // |t| must be greater than |current_time_|.  |planetarium_rotation| is the
  // value of KSP's |Planetarium.InverseRotAngle| at instant |t|, which provides
  // the rotation between the |World| axes and the |Barycentric| axes (we don't
  // use Planetarium.Rotation since it undergoes truncation to single-precision
  // even though it's a double-precision value).  Note that KSP's
  // |Planetarium.InverseRotAngle| is in degrees.
  virtual void AdvanceTime(Instant const& t, Angle const& planetarium_rotation);

  // Returns the degrees of freedom of the given part in |World|, assuming that
  // the origin of |World| is fixed at |bubble_barycentre_|.  The part must be
  // in a loaded vessel.
  virtual DegreesOfFreedom<World> GetPartActualDegreesOfFreedom(
      PartId part_id) const;

  virtual DegreesOfFreedom<Barycentric> GetBubbleBarycentre() const;

  // Returns the |World| degrees of freedom of the |Celestial| with the given
  // |Index|, identifying the origin of |World| with that of |Bubble|.
  virtual DegreesOfFreedom<World> CelestialWorldDegreesOfFreedom(
      Index const index) const;

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
      Index celestial_index) const;

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
      Index celestial_index,
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
  NewBarycentricRotatingNavigationFrame(Index primary_index,
                                        Index secondary_index) const;

  virtual not_null<std::unique_ptr<NavigationFrame>>
  NewBodyCentredBodyDirectionNavigationFrame(Index primary_index,
                                             Index secondary_index) const;

  virtual not_null<std::unique_ptr<NavigationFrame>>
  NewBodyCentredNonRotatingNavigationFrame(Index reference_body_index) const;

  virtual not_null<std::unique_ptr<NavigationFrame>>
  NewBodySurfaceNavigationFrame(Index reference_body_index) const;

  virtual void SetPlottingFrame(
      not_null<std::unique_ptr<NavigationFrame>> plotting_frame);
  virtual not_null<NavigationFrame const*> GetPlottingFrame() const;

  // The navball field at |current_time| for the current |plotting_frame_|.
  virtual std::unique_ptr<FrameField<World, Navball>> NavballFrameField(
      Position<World> const& sun_world_position) const;

  // The unit tangent, normal, or binormal vector to the trajectory of the
  // vessel with the given GUID in the current |plotting_frame_|.
  virtual Vector<double, World> VesselTangent(GUID const& vessel_guid) const;
  virtual Vector<double, World> VesselNormal(GUID const& vessel_guid) const;
  virtual Vector<double, World> VesselBinormal(GUID const& vessel_guid) const;

  virtual Velocity<World> VesselVelocity(GUID const& vessel_guid) const;

  // Coordinate transforms.
  virtual AffineMap<Barycentric, World, Length, OrthogonalMap>
  BarycentricToWorld(Position<World> const& sun_world_position) const;
  virtual OrthogonalMap<Barycentric, World> BarycentricToWorld() const;
  virtual OrthogonalMap<Barycentric, WorldSun> BarycentricToWorldSun() const;
  virtual AffineMap<World, Barycentric, Length, OrthogonalMap>
  WorldToBarycentric(Position<World> const& sun_world_position) const;
  virtual OrthogonalMap<World, Barycentric> WorldToBarycentric() const;

  virtual Instant GameEpoch() const;
  virtual bool MustRotateBodies() const;

  virtual Instant CurrentTime() const;

  // Must be called after initialization.
  virtual void WriteToMessage(not_null<serialization::Plugin*> message) const;
  static not_null<std::unique_ptr<Plugin>> ReadFromMessage(
      serialization::Plugin const& message);

 protected:
  // May be overriden in tests to inject a mock.
  virtual std::unique_ptr<Ephemeris<Barycentric>> NewEphemeris(
      std::vector<not_null<std::unique_ptr<MassiveBody const>>>&& bodies,
      std::vector<DegreesOfFreedom<Barycentric>> const& initial_state,
      Instant const& initial_time,
      Length const& fitting_tolerance,
      Ephemeris<Barycentric>::FixedStepParameters const& parameters);

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
         std::unique_ptr<Ephemeris<Barycentric>> ephemeris,
         Ephemeris<Barycentric>::FixedStepParameters const& history_parameters,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prolongation_parameters,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             prediction_parameters,
         Angle const& planetarium_rotation,
         Instant const& game_epoch,
         Instant const& current_time,
         Index sun_index,
         bool is_pre_cardano);

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

  Vector<double, World> FromVesselFrenetFrame(
      Vessel const& vessel,
      Vector<double, Frenet<Navigation>> const& vector) const;

  // Fill |celestials| using the |index| and |parent_index| fields found in
  // |celestial_messages| (which may be pre- or post-Bourbaki).
  template<typename T>
  static void ReadCelestialsFromMessages(
      Ephemeris<Barycentric> const& ephemeris,
      google::protobuf::RepeatedPtrField<T> const& celestial_messages,
      IndexToOwnedCelestial& celestials);

  // Computes a fingerprint for the parameters passed to
  // |InsertCelestialJacobiKeplerian|.
  static std::uint64_t FingerprintCelestialJacobiKeplerian(
      Index celestial_index,
      std::experimental::optional<Index> const& parent_index,
      std::experimental::optional<
          physics::KeplerianElements<Barycentric>> const& keplerian_elements,
      MassiveBody const& body);

  // Adds a part to a vessel, recording it in the appropriate map and setting up
  // a deletion callback.
  void AddPart(not_null<Vessel*> vessel,
               PartId part_id,
               std::string const& name,
               Mass mass,
               DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Whether |loaded_vessels_| contains |vessel|.
  bool is_loaded(not_null<Vessel*> vessel) const;
  // Whether |new_unloaded_vessels_| contains |vessel|.
  bool is_new_unloaded(not_null<Vessel*> vessel) const;

  GUIDToOwnedVessel vessels_;
  // For each part, the vessel that this part belongs to. The part is guaranteed
  // to be in the parts() map of the vessel, and owned by it.
  std::map<PartId, not_null<Vessel*>> part_id_to_vessel_;
  IndexToOwnedCelestial celestials_;

  // The vessels that will be kept during the next call to |AdvanceTime|.
  std::set<not_null<Vessel const*>> kept_vessels_;

  struct AbsoluteInitializationObjects final{
    IndexToMassiveBody bodies;
    IndexToDegreesOfFreedom initial_state;
  };
  std::experimental::optional<AbsoluteInitializationObjects>
      absolute_initialization_;

  struct HierarchicalInitializationObjects final {
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

  RotatingBody<Barycentric> const* main_body_ = nullptr;

  // Do not |erase| from this list, use |Vessel::clear_pile_up| instead.
  std::list<PileUp> pile_ups_;

  // The vessels that are currently loaded, i.e. in the physics bubble.
  std::set<not_null<Vessel*>> loaded_vessels_;
  // The vessels that were inserted unloaded and have yet to be collected into a
  // pile-up.
  std::set<not_null<Vessel*>> new_unloaded_vessels_;
  std::experimental::optional<DegreesOfFreedom<Barycentric>> bubble_barycentre_;

  // Compatibility.
  bool is_pre_cardano_ = false;

  friend class NavballFrameField;
  friend class TestablePlugin;
};

}  // namespace internal_plugin

using internal_plugin::Index;
using internal_plugin::Plugin;

}  // namespace ksp_plugin
}  // namespace principia
