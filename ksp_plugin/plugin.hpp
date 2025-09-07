#pragma once

#include <future>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "base/disjoint_sets.hpp"
#include "base/monostable.hpp"
#include "base/not_null.hpp"
#include "base/thread_pool.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/perspective.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "integrators/integrators.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/geometric_potential_plotter.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/planetarium.hpp"
#include "ksp_plugin/renderer.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/ephemeris.hpp"
#include "physics/frame_field.hpp"
#include "physics/hierarchical_system.hpp"
#include "physics/massive_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "physics/rotating_body.hpp"
#include "physics/tensors.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/astronomy.pb.h"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {

class TestablePlugin;

namespace _plugin {
namespace internal {

using namespace principia::base::_disjoint_sets;
using namespace principia::base::_monostable;
using namespace principia::base::_not_null;
using namespace principia::base::_thread_pool;
using namespace principia::geometry::_affine_map;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_perspective;
using namespace principia::geometry::_point;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::integrators::_integrators;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_geometric_potential_plotter;
using namespace principia::ksp_plugin::_identification;
using namespace principia::ksp_plugin::_pile_up;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::ksp_plugin::_renderer;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_body;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_segment;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_frame_field;
using namespace principia::physics::_hierarchical_system;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_tensors;
using namespace principia::physics::_trajectory;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// The index of a body in `FlightGlobals.Bodies`, obtained by
// `b.flightGlobalsIndex` in C#. We use this as a key in a map.
using Index = int;

class Plugin {
 public:
  Plugin() = delete;
  Plugin(Plugin const&) = delete;
  Plugin(Plugin&&) = delete;
  Plugin& operator=(Plugin const&) = delete;
  Plugin& operator=(Plugin&&) = delete;
  virtual ~Plugin();

  // Constructs a `Plugin`. The current time of that instance is
  // `solar_system_epoch`.  The angle between the axes of `World` and
  // `Barycentric` at `solar_system_epoch` is set to `planetarium_rotation`.
  // The epochs must be given in a format parseable by ParseTT.
  Plugin(std::string const& game_epoch,
         std::string const& solar_system_epoch,
         Angle const& planetarium_rotation);

  // Inserts a celestial body with index `celestial_index` and the given
  // `gravity_model` and `initial_state`.
  // If `parent_index` is null, inserts the sun, otherwise the parent of the new
  // body is the body with index `*parent_index`, which must already have been
  // inserted.
  // All the bodies must be inserted using the same method.
  virtual void InsertCelestialAbsoluteCartesian(
      Index celestial_index,
      std::optional<Index> const& parent_index,
      serialization::GravityModel::Body const& gravity_model,
      serialization::InitialState::Cartesian::Body const& initial_state);
  virtual void InsertCelestialJacobiKeplerian(
      Index celestial_index,
      std::optional<Index> const& parent_index,
      serialization::GravityModel::Body const& gravity_model,
      serialization::InitialState::Keplerian::Body const& initial_state);

  virtual void InitializeDownsamplingParameters(
      DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters const&
          downsampling_parameters);
  virtual void InitializeEphemerisParameters(
      Ephemeris<Barycentric>::AccuracyParameters const& accuracy_parameters,
      Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters);
  virtual void InitializeHistoryParameters(
      Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters);
  virtual void InitializePsychohistoryParameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const& parameters);
  // No setter for the default prediction parameters, as that default is always
  // overriden by the vessel-specific `SetPredictionAdaptiveStepParameters`.

  // Ends initialization.  The sun must have been inserted.
  virtual void EndInitialization();

  // And there shall, in that time, be rumors of things going astray, and there
  // will be a great confusion as to where things really are, and nobody will
  // really know where lieth those little things with the sort of raffia work
  // base, that has an attachment.
  virtual bool HasEncounteredApocalypse(std::string* details) const;

  // Sets the parent of the celestial body with index `celestial_index` to the
  // one with index `parent_index`. Both bodies must already have been
  // inserted. Must be called after initialization.
  // For a KSP `CelestialBody` `b`, the arguments correspond to
  // `b.flightGlobalsIndex`, `b.orbit.referenceBody.flightGlobalsIndex`.
  virtual void UpdateCelestialHierarchy(Index celestial_index,
                                        Index parent_index) const;

  // Sets the celestial whose axis of rotation will coincide with the `Alice`
  // z axis.
  virtual void SetMainBody(Index index);
  virtual Rotation<BodyWorld, World> CelestialRotation(Index index) const;
  virtual Rotation<CelestialSphere, World> CelestialSphereRotation() const;

  virtual Angle CelestialInitialRotation(Index celestial_index) const;
  virtual Time CelestialRotationPeriod(Index celestial_index) const;

  virtual void ClearWorldRotationalReferenceFrame();
  virtual void SetWorldRotationalReferenceFrame(Index celestial_index);

  virtual Index CelestialIndexOfBody(MassiveBody const& body) const;

  // Inserts a new vessel with GUID `vessel_guid` if it does not already exist,
  // and flags the vessel with GUID `vessel_guid` so it is kept when calling
  // `FreeVesselsAndPartsAndCollectPileUps`. The parent body for the vessel is
  // set to the one with index `parent_index`, which must have been inserted
  // during initialization.
  // Sets `inserted` to true if a new vessel was inserted, to false otherwise.
  // If `InsertOrKeepVessel` is called with `loaded=false`, and returns
  // `inserted=true`, `InsertUnloadedPart` must be called for its parts
  // before the call to `AdvanceTime`, giving the vessel an initial state.
  // For a KSP `Vessel` `v`, the arguments correspond to `v.id`,
  // `v.orbit.referenceBody.flightGlobalsIndex`, `v.loaded`.
  virtual void InsertOrKeepVessel(GUID const& vessel_guid,
                                  std::string const& vessel_name,
                                  Index parent_index,
                                  bool loaded,
                                  bool& inserted);

  // Adds a part with the given `part_id` to the vessel with the given `GUID`,
  // which must be unloaded, putting the part at the given offset from the
  // parent body of the vessel.  The part is given unit mass; this does not
  // matter, since the `PileUp` will be deformed when it first loads anyway.
  virtual void InsertUnloadedPart(
      PartId part_id,
      std::string const& name,
      GUID const& vessel_guid,
      RelativeDegreesOfFreedom<AliceSun> const& from_parent);

  // Inserts a new part with the given ID if it does not already exist, and
  // flags the part with the given ID so it is kept when calling
  // `FreeVesselsAndPartsAndCollectPileUps`.
  // The part is created in the given `vessel`, or if it already existed in
  // another vessel, is moved to that one.  If the part is created, its degrees
  // of freedom are set using those given and the `main_body_index`; otherwise
  // these three parameters are ignored.
  virtual void InsertOrKeepLoadedPart(
      PartId part_id,
      std::string const& name,
      Mass const& mass,
      Position<EccentricPart> const& centre_of_mass,
      InertiaTensor<RigidPart> const& inertia_tensor,
      bool is_solid_rocket_motor,
      GUID const& vessel_guid,
      Index main_body_index,
      DegreesOfFreedom<World> const& main_body_degrees_of_freedom,
      RigidMotion<EccentricPart, World> const& part_rigid_motion,
      Time const& Δt);

  // Calls `apply_intrinsic_force` and `apply_intrinsic_torque` on the
  // relevant part, which must be in a loaded vessel.
  virtual void ApplyPartIntrinsicForce(
      PartId part_id,
      Vector<Force, World> const& force) const;
  virtual void ApplyPartIntrinsicForceAtPosition(
      PartId part_id,
      Vector<Force, World> const& force,
      Displacement<World> const& lever_arm) const;
  virtual void ApplyPartIntrinsicTorque(
      PartId part_id,
      Bivector<Torque, World> const& torque) const;

  virtual bool PartIsTruthful(PartId part_id) const;

  // Calls `MakeSingleton` for all parts in loaded vessels, enabling the use of
  // union-find for pile up construction.  This must be called after the calls
  // to `ApplyPartIntrinsicForce`, and before the calls to
  // `ReportGroundCollision` or `ReportPartCollision`.
  virtual void PrepareToReportCollisions();

  // Notifies `this` that the given part is touching the ground.
  virtual void ReportGroundCollision(PartId part) const;

  // Notifies `this` that the given parts are touching, and should gravitate
  // as part of a single rigid body.
  virtual void ReportPartCollision(PartId part1, PartId part2) const;

  // Destroys the vessels for which `InsertOrKeepVessel` has not been called
  // since the last call to `FreeVesselsAndCollectPileUps`, as well as the
  // vessels which transitively touch the ground.  Destroys the parts in loaded
  // vessels for which `InsertOrKeepLoadedPart` has not been called.  Updates
  // the list of `pile_ups_` according to the reported collisions.
  virtual void FreeVesselsAndPartsAndCollectPileUps(Time const& Δt);

  // Calls `SetPartApparentRigidMotion` on the pile-up containing the relevant
  // part.  This part must be in a loaded vessel.
  virtual void SetPartApparentRigidMotion(
      PartId part_id,
      RigidMotion<EccentricPart, ApparentWorld> const& rigid_motion);

  // Returns the motion of the given part in `World`, assuming that
  // the origin of `World` is fixed at the centre of mass of the
  // `part_at_origin`.
  virtual RigidMotion<EccentricPart, World> GetPartActualMotion(
      PartId part_id,
      RigidMotion<Barycentric, World> const& barycentric_to_world) const;

  // Returns the `World` degrees of freedom of the `Celestial` with the given
  // `Index`, identifying the origin of `World` with the centre of mass of the
  // `Part` with the given `PartId`.
  virtual DegreesOfFreedom<World> CelestialWorldDegreesOfFreedom(
      Index index,
      RigidMotion<Barycentric, World> const& barycentric_to_world,
      Instant const& time) const;

  virtual RigidMotion<Barycentric, World> BarycentricToWorld(
      bool reference_part_is_unmoving,
      PartId reference_part_id,
      std::optional<Position<World>> const& main_body_centre) const;

  // Simulates the system until instant `t`.  Sets `current_time_` to `t`.
  // Must be called after initialization.
  // Clears the intrinsic force on all loaded parts.
  // `t` must be greater than `current_time_`.  `planetarium_rotation` is the
  // value of KSP's `Planetarium.InverseRotAngle` at instant `t`, which provides
  // the rotation between the `World` axes and the `Barycentric` axes (we don't
  // use Planetarium.Rotation since it undergoes truncation to single-precision
  // even though it's a double-precision value).  Note that KSP's
  // `Planetarium.InverseRotAngle` is in degrees.
  virtual void AdvanceTime(Instant const& t, Angle const& planetarium_rotation);

  // Advances time to `current_time_` for all pile ups that are not already
  // there, filling the tails of all their parts up to that instant; then
  // advances time on all vessels that are not yet at `current_time_`.  Inserts
  // the set of vessels that have collided with a celestial into
  // `collided_vessels`.
  virtual void CatchUpLaggingVessels(VesselSet& collided_vessels);

  // Advances time to `current_time_` on the pile up containing the given
  // vessel if the pile up is not there already, and advances time to
  // `current_time_` on that vessel.  This operation is asynchronous: the caller
  // must wait on the returned future before using the trajectories of the
  // vessel.  The caller must ensure that the vessels don't change while this
  // method is running.
  virtual not_null<std::unique_ptr<PileUpFuture>> CatchUpVessel(
      GUID const& vessel_guid);

  // Waits for the `future` to return and inserts the set of vessels that have
  // collided with a celestial into `collided_vessels`.
  virtual void WaitForVesselToCatchUp(PileUpFuture& pile_up_future,
                                      VesselSet& collided_vessels);

  // Returns the displacement and velocity of the vessel with GUID `vessel_guid`
  // relative to its parent at current time. For a KSP `Vessel` `v`, the
  // argument corresponds to  `v.id.ToString()`, the return value to
  // `{v.orbit.pos, v.orbit.vel}`.
  // A vessel with GUID `vessel_guid` must have been inserted and kept. Must
  // be called after initialization.
  virtual RelativeDegreesOfFreedom<AliceSun> VesselFromParent(
      Index parent_index,
      GUID const& vessel_guid) const;

  // Returns the displacement and velocity of the celestial at index
  // `celestial_index` relative to its parent at current time. For a KSP
  // `CelestialBody` `b`, the argument corresponds to `b.flightGlobalsIndex`,
  // the return value to `{b.orbit.pos, b.orbit.vel}`.
  // A celestial with index `celestial_index` must have been inserted, and it
  // must not be the sun. Must be called after initialization.
  virtual RelativeDegreesOfFreedom<AliceSun> CelestialFromParent(
      Index celestial_index) const;

  virtual void SetPredictionAdaptiveStepParameters(
      GUID const& vessel_guid,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          prediction_adaptive_step_parameters) const;

  // Updates the prediction for the vessels with guids in `vessel_guids`.
  void UpdatePrediction(std::vector<GUID> const& vessel_guids) const;

  virtual void CreateFlightPlan(GUID const& vessel_guid,
                                Instant const& final_time,
                                Mass const& initial_mass) const;

  // If there is a target vessel, increases the size of the prediction to try
  // and cover the actual final time of the flight plan.
  virtual void ExtendPredictionForFlightPlan(GUID const& vessel_guid) const;

  // Computes the apsides of the trajectory defined by `begin` and `end` with
  // respect to the celestial with index `celestial_index`.
  virtual void ComputeAndRenderApsides(
      Index celestial_index,
      Trajectory<Barycentric> const& trajectory,
      DiscreteTrajectory<Barycentric>::iterator const& begin,
      DiscreteTrajectory<Barycentric>::iterator const& end,
      Instant const& t_max,
      Position<World> const& sun_world_position,
      int max_points,
      DiscreteTrajectory<World>& apoapsides,
      DiscreteTrajectory<World>& periapsides) const;

  // Computes the first collision between the trajectory defined by `begin` and
  // `end` and the celestial with index `celestial_index`.
  virtual std::optional<DiscreteTrajectory<World>::value_type>
  ComputeAndRenderFirstCollision(
      Index celestial_index,
      Trajectory<Barycentric> const& trajectory,
      DiscreteTrajectory<Barycentric>::iterator const& begin,
      DiscreteTrajectory<Barycentric>::iterator const& end,
      Position<World> const& sun_world_position,
      int max_points,
      std::function<Length(Angle const& latitude,
                           Angle const& longitude)> const& radius) const;

  // Computes the closest approaches of the trajectory defined by `begin` and
  // `end` with respect to the trajectory of the targetted vessel.
  virtual void ComputeAndRenderClosestApproaches(
      Trajectory<Barycentric> const& trajectory,
      DiscreteTrajectory<Barycentric>::iterator const& begin,
      DiscreteTrajectory<Barycentric>::iterator const& end,
      Position<World> const& sun_world_position,
      int max_points,
      DiscreteTrajectory<World>& closest_approaches) const;

  // Computes the nodes of the trajectory defined by `begin` and `end` with
  // respect to plane of the trajectory of the targetted vessel.
  virtual void ComputeAndRenderNodes(
      DiscreteTrajectory<Barycentric>::iterator const& begin,
      DiscreteTrajectory<Barycentric>::iterator const& end,
      Instant const& t_max,
      Position<World> const& sun_world_position,
      int max_points,
      std::vector<Renderer::Node>& ascending,
      std::vector<Renderer::Node>& descending) const;

  virtual bool HasCelestial(Index index) const;
  virtual Celestial const& GetCelestial(Index index) const;

  virtual bool HasVessel(GUID const& vessel_guid) const;
  virtual not_null<Vessel*> GetVessel(GUID const& vessel_guid) const;

  virtual void ClearOrbitAnalysersOfVesselsOtherThan(Vessel const& vessel);

  virtual not_null<std::unique_ptr<Planetarium>> NewPlanetarium(
      Planetarium::Parameters const& parameters,
      Perspective<Navigation, Camera> const& perspective,
      Planetarium::PlottingToScaledSpaceConversion plotting_to_scaled_space)
      const;

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

  virtual not_null<std::unique_ptr<PlottingFrame>>
  NewRotatingPulsatingPlottingFrame(
      std::vector<Index> const& primary_indices,
      std::vector<Index> const& secondary_indices) const;

  virtual void SetTargetVessel(GUID const& vessel_guid,
                               Index reference_body_index);

  // The navball field at `current_time` for the current `plotting_frame_`.
  virtual std::unique_ptr<FrameField<World, Navball>> NavballFrameField(
      Position<World> const& sun_world_position) const;

  // The unit tangent, normal, or binormal vector to the trajectory of the
  // vessel with the given GUID in the current `plotting_frame_`.
  virtual Vector<double, World> VesselTangent(GUID const& vessel_guid) const;
  virtual Vector<double, World> VesselNormal(GUID const& vessel_guid) const;
  virtual Vector<double, World> VesselBinormal(GUID const& vessel_guid) const;

  // TODO(egg): UnmanageableVesselTangent, Normal, Binormal.

  // Takes degrees of freedom relative to the celestial with the given index,
  // and returns the velocity in the plotting frame expressed in the coordinates
  // of `World`.  This is used to display the velocity of a vessel not known to
  // the plugin.
  virtual Velocity<World> UnmanageableVesselVelocity(
      RelativeDegreesOfFreedom<AliceSun> const& degrees_of_freedom,
      Index parent_index) const;
  // Same as `UnmanageableVesselVelocity`, but uses the known degrees of freedom
  // of a vessel in `vessels_`.
  virtual Velocity<World> VesselVelocity(GUID const& vessel_guid) const;

  void RequestReanimation(Instant const& desired_t_min) const;

  virtual Instant GameEpoch() const;

  virtual Instant CurrentTime() const;

  // The rotation between the `AliceWorld` basis at `current_time_` and the
  // `Barycentric` axes. Since `AliceSun` is not a rotating reference frame,
  // this change of basis is all that's required to convert relative velocities
  // or displacements between simultaneous events.
  virtual Rotation<Barycentric, AliceSun> const& PlanetariumRotation() const;
  virtual Rotation<CameraCompensatedReference, CameraReference> const&
  CameraCompensation() const;

  virtual Renderer& renderer();
  virtual Renderer const& renderer() const;

  virtual GeometricPotentialPlotter& geometric_potential_plotter();
  virtual GeometricPotentialPlotter const& geometric_potential_plotter() const;

  // Must be called after initialization.
  virtual void WriteToMessage(not_null<serialization::Plugin*> message) const;
  static not_null<std::unique_ptr<Plugin>> ReadFromMessage(
      serialization::Plugin const& message);

 private:
  using GUIDToOwnedVessel = std::map<GUID, not_null<std::unique_ptr<Vessel>>>;
  using IndexToOwnedCelestial =
      std::map<Index, not_null<std::unique_ptr<Celestial>>>;
  using NewtonianMotionEquation =
      Ephemeris<Barycentric>::NewtonianMotionEquation;

  // This constructor should only be used during deserialization.  The plugin is
  // not fully initialized when it returns.
  Plugin(Ephemeris<Barycentric>::FixedStepParameters history_parameters,
         Ephemeris<Barycentric>::AdaptiveStepParameters
             psychohistory_parameters);

  void InitializeIndices(
      std::string const& name,
      Index celestial_index,
      std::optional<Index> const& parent_index);

  // Computes the value returned by `PlanetariumRotation`.  Must be called
  // whenever `main_body_` or `planetarium_rotation_` changes.
  void UpdatePlanetariumRotation();

  Velocity<World> VesselVelocity(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom) const;

  // Fill `celestials` using the `index` and `parent_index` fields found in
  // `celestial_messages`.
  template<typename T>
  static void ReadCelestialsFromMessages(
      Ephemeris<Barycentric> const& ephemeris,
      google::protobuf::RepeatedPtrField<T> const& celestial_messages,
      IndexToOwnedCelestial& celestials,
      std::map<std::string, Index>& name_to_index);

  // Constructs a part using the constructor arguments, and add it to a vessel,
  // recording it in the appropriate map and setting up a deletion callback.
  template<typename... Args>
  void AddPart(not_null<Vessel*> vessel,
               PartId part_id,
               std::string const& name,
               Args... args);

  // Whether `loaded_vessels_` contains `vessel`.
  bool is_loaded(not_null<Vessel*> vessel) const;

  // Initialization objects.
  Monostable initializing_;
  bool uses_correct_sin_cos_ = true;
  serialization::GravityModel gravity_model_;
  serialization::InitialState initial_state_;
  std::map<std::string, Index> name_to_index_;
  std::map<Index, std::string> index_to_name_;
  std::uint64_t system_fingerprint_ = 0;
  std::map<Index, std::optional<Index>> parents_;
  // The ephemeris is only constructed once, so this is an initialization
  // object.  The other parameters must be persisted to create new vessels.
  // Since this is not persisted directly, it is optional so that it can be null
  // in a deserialized object.
  std::optional<Ephemeris<Barycentric>::AccuracyParameters>
      ephemeris_accuracy_parameters_;
  std::optional<Ephemeris<Barycentric>::FixedStepParameters>
      ephemeris_fixed_step_parameters_;

  GUIDToOwnedVessel vessels_;
  // For each part, the vessel that this part belongs to. The part is guaranteed
  // to be in the parts() map of the vessel, and owned by it.
  absl::flat_hash_map<PartId, not_null<Vessel*>> part_id_to_vessel_;
  IndexToOwnedCelestial celestials_;

  // Not null after initialization.
  std::unique_ptr<Ephemeris<Barycentric>> ephemeris_;

  // The parameters for computing the various trajectories.
  DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
      history_downsampling_parameters_;
  Ephemeris<Barycentric>::FixedStepParameters history_fixed_step_parameters_;
  Ephemeris<Barycentric>::AdaptiveStepParameters psychohistory_parameters_;

  // The thread pool for advancing vessels.
  ThreadPool<absl::Status> vessel_thread_pool_;

  Angle planetarium_rotation_;
  std::optional<Rotation<Barycentric, AliceSun>> cached_planetarium_rotation_;
  std::optional<Rotation<CameraCompensatedReference, CameraReference>>
      camera_compensation_;
  // The game epoch in real time.
  Instant game_epoch_;
  // The current in-game universal time.
  Instant current_time_;

  Celestial* sun_ = nullptr;  // Not owning, not null after initialization.

  // Not null after initialization.
  std::unique_ptr<Renderer> renderer_;

  RotatingBody<Barycentric> const* main_body_ = nullptr;
  AngularVelocity<Barycentric> angular_velocity_of_world_;

  // Do not `erase` from this list, use `Part::reset_containing_pile_up` instead
  // and the pile-up will remove itself once no part owns it.  The elements are
  // not `not_null<>` because we temporarily need to insert null pointers.
  std::list<PileUp*> pile_ups_;

  // The vessels that are currently loaded, i.e. in the physics bubble.
  VesselSet loaded_vessels_;
  // The vessels that will be kept during the next call to `AdvanceTime`.
  VesselConstSet kept_vessels_;
  // Contains the adaptive step parameters for the vessel that existed in the
  // past but are no longer known to the plugin.  Useful to avoid losing the
  // parameters, e.g., when a vessel hits the ground.
  // NOTE(phl): This is a leaky map, in the sense that we don't remove deleted
  // vessels from it.  Hopefully it's small enough that we don't care.
  std::map<GUID, Ephemeris<Barycentric>::AdaptiveStepParameters>
  zombie_prediction_adaptive_step_parameters_;

  std::optional<GeometricPotentialPlotter> geometric_potential_plotter_;

  friend class NavballFrameField;
  friend class ksp_plugin::TestablePlugin;
};

}  // namespace internal

using internal::Index;
using internal::Plugin;

}  // namespace _plugin
}  // namespace ksp_plugin
}  // namespace principia
