
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <set>

#include "base/map_util.hpp"
#include "base/not_null.hpp"
#include "base/optional_logging.hpp"
#include "base/unique_ptr_logging.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/identity.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/permutation.hpp"
#include "glog/logging.h"
#include "glog/stl_logging.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/barycentric_rotating_dynamic_frame_body.hpp"
#include "physics/body_centered_non_rotating_dynamic_frame.hpp"
#include "physics/dynamic_frame.hpp"

namespace principia {
namespace ksp_plugin {

using base::FindOrDie;
using base::make_not_null_unique;
using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Identity;
using geometry::Normalize;
using geometry::Permutation;
using geometry::Sign;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using physics::BarycentricRotatingDynamicFrame;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DynamicFrame;
using physics::Frenet;
using quantities::Force;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using ::operator<<;

namespace {

// The map between the vector spaces of |World| and |AliceWorld|.
Permutation<World, AliceWorld> const kWorldLookingGlass(
    Permutation<World, AliceWorld>::CoordinatePermutation::XZY);

// The map between the vector spaces of |WorldSun| and |AliceSun|.
Permutation<WorldSun, AliceSun> const kSunLookingGlass(
    Permutation<WorldSun, AliceSun>::CoordinatePermutation::XZY);

Time const kStep = 45 * Minute;
Length const kFittingTolerance = 1 * Milli(Metre);

}  // namespace

Plugin::Plugin(Instant const& initial_time,
               Angle const& planetarium_rotation)
    : bubble_(make_not_null_unique<PhysicsBubble>()),
      history_integrator_(
          McLachlanAtela1992Order5Optimal<Position<Barycentric>>()),
      prolongation_integrator_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>()),
      prediction_integrator_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>()),
      planetarium_rotation_(planetarium_rotation),
      current_time_(initial_time),
      history_time_(initial_time) {}

// The three celestial-inserting functions log at the INFO level.  This is not
// a concern since they are only called at initialization, which should happen
// once per game.

void Plugin::InsertSun(Index const celestial_index,
                       GravitationalParameter const& gravitational_parameter) {
  CHECK(initializing_) << "Celestial bodies should be inserted before the end "
                       << "of initialization";
  CHECK(!absolute_initialization_);
  CHECK(!hierarchical_initialization_);
  LOG(INFO) << __FUNCTION__ << "\n"
            << NAMED(celestial_index) << "\n"
            << NAMED(gravitational_parameter);
  auto sun = make_not_null_unique<MassiveBody>(gravitational_parameter);
  auto const unowned_sun = sun.get();
  hierarchical_initialization_.emplace(std::move(sun));
  hierarchical_initialization_->indices_to_bodies[celestial_index] =
      unowned_sun;
  hierarchical_initialization_->parents[celestial_index] =
      std::experimental::nullopt;
}

void Plugin::InsertCelestialAbsoluteCartesian(
    Index const celestial_index,
    std::experimental::optional<Index> const& parent_index,
    DegreesOfFreedom<Barycentric> const& initial_state,
    base::not_null<std::unique_ptr<MassiveBody const>> body) {
  CHECK(initializing_) << "Celestial bodies should be inserted before the end "
                       << "of initialization";
  CHECK(!hierarchical_initialization_);
  LOG(INFO) << __FUNCTION__ << "\n"
            << NAMED(celestial_index) << "\n"
            << NAMED(parent_index) << "\n"
            << NAMED(initial_state) << "\n"
            << NAMED(body);
  if (!absolute_initialization_) {
    absolute_initialization_.emplace();
  }
  auto const inserted = celestials_.emplace(
      celestial_index,
      std::make_unique<Celestial>(body.get()));
  CHECK(inserted.second) << "Body already exists at index " << celestial_index;
  not_null<Celestial*> const celestial = inserted.first->second.get();
  absolute_initialization_->bodies.emplace(celestial_index, std::move(body));
  if (parent_index) {
    not_null<Celestial const*> parent =
        FindOrDie(celestials_, *parent_index).get();
    celestial->set_parent(parent);
  } else {
    CHECK(sun_ == nullptr);
    sun_ = celestial;
  }
  absolute_initialization_->initial_state.emplace(celestial_index,
                                                  initial_state);
}

void Plugin::InsertCelestialJacobiKeplerian(
    Index const celestial_index,
    Index const parent_index,
    KeplerianElements<Barycentric> const& keplerian_elements,
    base::not_null<std::unique_ptr<MassiveBody>> body) {
  LOG(INFO) << __FUNCTION__ << "\n"
            << NAMED(celestial_index) << "\n"
            << NAMED(parent_index) << "\n"
            << NAMED(keplerian_elements) << "\n"
            << NAMED(body);
  CHECK(initializing_);
  CHECK(hierarchical_initialization_);
  bool inserted =
      hierarchical_initialization_->parents.emplace(celestial_index,
                                                    parent_index).second;
  inserted &=
      hierarchical_initialization_->
          indices_to_bodies.emplace(celestial_index, body.get()).second;
  CHECK(inserted);
  hierarchical_initialization_->system.Add(
      std::move(body),
      hierarchical_initialization_->indices_to_bodies[parent_index],
      keplerian_elements);
}

void Plugin::EndInitialization() {
  CHECK(initializing_);
  if (hierarchical_initialization_) {
    HierarchicalSystem<Barycentric>::BarycentricSystem system =
        hierarchical_initialization_->system.ConsumeBarycentricSystem();
    std::map<not_null<MassiveBody const*>, Index> bodies_to_indices;
    for (auto const& index_body :
         hierarchical_initialization_->indices_to_bodies) {
      bodies_to_indices[index_body.second] = index_body.first;
    }
    auto const parents = std::move(hierarchical_initialization_->parents);
    hierarchical_initialization_ = std::experimental::nullopt;
    for (int i = 0; i < system.bodies.size(); ++i) {
      Index const celestial_index = bodies_to_indices[system.bodies[i].get()];
      InsertCelestialAbsoluteCartesian(
          celestial_index,
          FindOrDie(parents, celestial_index),
          system.degrees_of_freedom[i],
          std::move(system.bodies[i]));
    }
  }
  CHECK(absolute_initialization_);
  CHECK_NOTNULL(sun_);
  initializing_.Flop();

  InitializeEphemerisAndSetCelestialTrajectories();
}

void Plugin::UpdateCelestialHierarchy(Index const celestial_index,
                                      Index const parent_index) const {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(celestial_index) << '\n' << NAMED(parent_index);
  CHECK(!initializing_);
  FindOrDie(celestials_, celestial_index)->set_parent(
      FindOrDie(celestials_, parent_index).get());
}

bool Plugin::InsertOrKeepVessel(GUID const& vessel_guid,
                                Index const parent_index) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(vessel_guid) << '\n' << NAMED(parent_index);
  CHECK(!initializing_);
  not_null<Celestial const*> parent =
      FindOrDie(celestials_, parent_index).get();
  auto inserted = vessels_.emplace(vessel_guid,
                                   make_not_null_unique<Vessel>(parent));
  not_null<Vessel*> const vessel = inserted.first->second.get();
  kept_vessels_.emplace(vessel);
  vessel->set_parent(parent);
  LOG_IF(INFO, inserted.second) << "Inserted vessel with GUID " << vessel_guid
                                << " at " << vessel;
  VLOG(1) << "Parent of vessel with GUID " << vessel_guid <<" is at index "
          << parent_index;
  return inserted.second;
}

void Plugin::SetVesselStateOffset(
    GUID const& vessel_guid,
    RelativeDegreesOfFreedom<AliceSun> const& from_parent) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(vessel_guid) << '\n' << NAMED(from_parent);
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(!vessel->is_initialized())
      << "Vessel with GUID " << vessel_guid << " already has a trajectory";
  LOG(INFO) << "Initial |{orbit.pos, orbit.vel}| for vessel with GUID "
            << vessel_guid << ": " << from_parent;
  RelativeDegreesOfFreedom<Barycentric> const relative =
      PlanetariumRotation().Inverse()(from_parent);
  LOG(INFO) << "In barycentric coordinates: " << relative;
  ephemeris_->Prolong(current_time_);
  vessel->CreateProlongation(
      current_time_,
      vessel->parent()->current_degrees_of_freedom(current_time_) + relative);
  auto const inserted = unsynchronized_vessels_.emplace(vessel.get());
  CHECK(inserted.second);
}

void Plugin::AdvanceTime(Instant const& t, Angle const& planetarium_rotation) {
  VLOG(1) << __FUNCTION__ << '\n'
          << NAMED(t) << '\n' << NAMED(planetarium_rotation);
  CHECK(!initializing_);
  CHECK_GT(t, current_time_);
  CleanUpVessels();
  ephemeris_->Prolong(t);
  bubble_->Prepare(BarycentricToWorldSun(), current_time_, t);
  Trajectories synchronized_histories =
      SynchronizedHistories();
  bool advanced_history_time = false;
  if (synchronized_histories.empty()) {
    // Synchronize everything now, nothing is currently synchronized.
    history_time_ = t;
    advanced_history_time = true;
  } else if (history_time_ + Δt_ <= t) {
    // The histories are far enough behind that we can advance them at least one
    // step and reset the prolongations.
    EvolveHistories(t, synchronized_histories);
    advanced_history_time = true;
  }
  if (advanced_history_time) {
    // TODO(egg): I think |!bubble_->empty()| => |has_dirty_vessels()|.
    if (has_unsynchronized_vessels() ||
        has_dirty_vessels() ||
        !bubble_->empty()) {
      SynchronizeNewVesselsAndCleanDirtyVessels();
    }
    ResetProlongations();
  }
  if (history_time_ < t) {
    EvolveProlongationsAndBubble(t);
  }
  VLOG(1) << "Time has been advanced" << '\n'
          << "from : " << current_time_ << '\n'
          << "to   : " << t;
  current_time_ = t;
  planetarium_rotation_ = planetarium_rotation;
}

void Plugin::ForgetAllHistoriesBefore(Instant const& t) const {
  CHECK(!initializing_);
  CHECK_LT(t, history_time_);

  // If |t| is after the start of the flight plan of some vessel, adjust it
  // to avoid destroying the trajectories of that flight plan.
  Instant forgetting_time = t;
  for (auto const& pair : vessels_) {
    not_null<std::unique_ptr<Vessel>> const& vessel = pair.second;
    if (vessel->has_flight_plan()) {
      forgetting_time =
          std::min(forgetting_time,
                   vessel->flight_plan()->initial_time());
    }
  }

  ephemeris_->ForgetBefore(forgetting_time);
  for (auto const& pair : vessels_) {
    not_null<std::unique_ptr<Vessel>> const& vessel = pair.second;
    // Only forget the synchronized vessels, the others don't have an history.
    if (unsynchronized_vessels_.count(vessel.get()) == 0) {
      vessel->mutable_history()->ForgetBefore(forgetting_time);
    }
  }
}

RelativeDegreesOfFreedom<AliceSun> Plugin::VesselFromParent(
    GUID const& vessel_guid) const {
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << vessel_guid
                                  << " was not given an initial state";
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      vessel->prolongation().last().degrees_of_freedom() -
      vessel->parent()->current_degrees_of_freedom(current_time_);
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  VLOG(1) << "Vessel with GUID " << vessel_guid
          << " is at parent degrees of freedom + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

RelativeDegreesOfFreedom<AliceSun> Plugin::CelestialFromParent(
    Index const celestial_index) const {
  CHECK(!initializing_);
  ephemeris_->Prolong(current_time_);
  Celestial const& celestial = *FindOrDie(celestials_, celestial_index);
  CHECK(celestial.has_parent())
      << "Body at index " << celestial_index << " is the sun";
  RelativeDegreesOfFreedom<Barycentric> const barycentric_result =
      celestial.current_degrees_of_freedom(CurrentTime()) -
      celestial.parent()->current_degrees_of_freedom(CurrentTime());
  RelativeDegreesOfFreedom<AliceSun> const result =
      PlanetariumRotation()(barycentric_result);
  VLOG(1) << "Celestial at index " << celestial_index
          << " is at parent degrees of freedom + " << barycentric_result
          << " Barycentre (" << result << " AliceSun)";
  return result;
}

void Plugin::UpdatePrediction(GUID const& vessel_guid) const {
  CHECK(!initializing_);
  find_vessel_by_guid_or_die(vessel_guid)->UpdatePrediction(
      ephemeris_.get(),
      prediction_integrator_,
      current_time_ + prediction_length_,
      prediction_length_tolerance_,
      prediction_speed_tolerance_);
}

void Plugin::CreateFlightPlan(GUID const& vessel_guid,
                              Instant const& final_time,
                              Mass const& initial_mass) const {
  CHECK(!initializing_);
  find_vessel_by_guid_or_die(vessel_guid)->CreateFlightPlan(
      final_time,
      initial_mass,
      ephemeris_.get(),
      prediction_integrator_,
      prediction_length_tolerance_,
      prediction_speed_tolerance_);
}

RenderedTrajectory<World> Plugin::RenderedVesselTrajectory(
    GUID const& vessel_guid,
    Position<World> const& sun_world_position) const {
  CHECK(!initializing_);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK(vessel->is_initialized());
  VLOG(1) << "Rendering a trajectory for the vessel with GUID " << vessel_guid;
  if (!vessel->is_synchronized()) {
    // TODO(egg): We render neither unsynchronized histories nor prolongations
    // at the moment.
    VLOG(1) << "Returning an empty trajectory";
    return RenderedTrajectory<World>();
  }

  return RenderedTrajectoryFromIterators(vessel->history().Begin(),
                                         vessel->history().End(),
                                         sun_world_position);
}

bool Plugin::HasPrediction(GUID const& vessel_guid) const {
  return find_vessel_by_guid_or_die(vessel_guid)->has_prediction();
}

RenderedTrajectory<World> Plugin::RenderedPrediction(
    GUID const& vessel_guid,
    Position<World> const& sun_world_position) const {
  CHECK(!initializing_);
  Vessel const& vessel = *find_vessel_by_guid_or_die(vessel_guid);
  RenderedTrajectory<World> result =
      RenderedTrajectoryFromIterators(vessel.prediction().Fork(),
                                      vessel.prediction().End(),
                                      sun_world_position);
  return result;
}

RenderedTrajectory<World> Plugin::RenderedTrajectoryFromIterators(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position) const {
  RenderedTrajectory<World> result;
  auto const to_world =
      AffineMap<Barycentric, World, Length, OrthogonalMap>(
          sun_->current_position(current_time_),
          sun_world_position,
          OrthogonalMap<WorldSun, World>::Identity() * BarycentricToWorldSun());

  // Compute the trajectory in the navigation frame.
  DiscreteTrajectory<Navigation> intermediate_trajectory;
  for (auto it = begin; it != end; ++it) {
    intermediate_trajectory.Append(
        it.time(),
        plotting_frame_->ToThisFrameAtTime(it.time())(
            it.degrees_of_freedom()));
  }

  // Render the trajectory at current time in |World|.
  DiscreteTrajectory<Navigation>::Iterator initial_it =
      intermediate_trajectory.Begin();
  DiscreteTrajectory<Navigation>::Iterator const intermediate_end =
      intermediate_trajectory.End();
  auto from_navigation_frame_to_world_at_current_time =
      to_world *
          plotting_frame_->
              FromThisFrameAtTime(current_time_).rigid_transformation();
  if (initial_it != intermediate_end) {
    for (auto final_it = initial_it;
         ++final_it, final_it != intermediate_end;
         initial_it = final_it) {
      result.emplace_back(from_navigation_frame_to_world_at_current_time(
                              initial_it.degrees_of_freedom().position()),
                          from_navigation_frame_to_world_at_current_time(
                              final_it.degrees_of_freedom().position()));
    }
  }
  VLOG(1) << "Returning a " << result.size() << "-segment trajectory";
  return result;
}

void Plugin::SetPredictionLength(Time const& t) {
  prediction_length_ = t;
}

void Plugin::SetPredictionLengthTolerance(Length const& l) {
  prediction_length_tolerance_ = l;
}

void Plugin::SetPredictionSpeedTolerance(Speed const& v) {
  prediction_speed_tolerance_ = v;
}

bool Plugin::HasVessel(GUID const& vessel_guid) const {
  return vessels_.find(vessel_guid) != vessels_.end();
}

not_null<Vessel*> Plugin::GetVessel(GUID const& vessel_guid) const {
  CHECK(!initializing_);
  return find_vessel_by_guid_or_die(vessel_guid).get();
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBodyCentredNonRotatingNavigationFrame(
    Index const reference_body_index) const {
  CHECK(!initializing_);
  Celestial const& reference_body =
      *FindOrDie(celestials_, reference_body_index);
  return make_not_null_unique<
      BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          reference_body.body());
}

not_null<std::unique_ptr<NavigationFrame>>
Plugin::NewBarycentricRotatingNavigationFrame(
    Index const primary_index,
    Index const secondary_index) const {
  CHECK(!initializing_);
  // TODO(egg): these should be const, use a custom comparator in the map.
  Celestial const& primary = *FindOrDie(celestials_, primary_index);
  Celestial const& secondary = *FindOrDie(celestials_, secondary_index);
  return make_not_null_unique<
      BarycentricRotatingDynamicFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          primary.body(),
          secondary.body());
}

void Plugin::SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) {
  plotting_frame_ = std::move(plotting_frame);
}

not_null<NavigationFrame const*> Plugin::GetPlottingFrame() const {
  return plotting_frame_.get();
}

Position<World> Plugin::PlotBarycentricPosition(
    Instant const& t,
    Position<Barycentric> const& position,
    Position<World> const& sun_world_position) const {
  auto const barycentric_to_world =
      AffineMap<Barycentric, World, Length, OrthogonalMap>(
          sun_->current_position(current_time_),
          sun_world_position,
          OrthogonalMap<WorldSun, World>::Identity() * BarycentricToWorldSun());
  auto const barycentric_to_navigation_at_t =
      plotting_frame_->ToThisFrameAtTime(t).rigid_transformation();
  auto const navigation_frame_to_world_at_current_time =
      barycentric_to_world *
          plotting_frame_->
              FromThisFrameAtTime(current_time_).rigid_transformation();
  return navigation_frame_to_world_at_current_time(
             barycentric_to_navigation_at_t(position));
}

void Plugin::AddVesselToNextPhysicsBubble(
    GUID const& vessel_guid,
    std::vector<IdAndOwnedPart> parts) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel_guid) << '\n' << NAMED(parts);
  not_null<std::unique_ptr<Vessel>> const& vessel =
      find_vessel_by_guid_or_die(vessel_guid);
  CHECK_LT(0, kept_vessels_.count(vessel.get()));
  dirty_vessels_.insert(vessel.get());
  bubble_->AddVesselToNext(vessel.get(), std::move(parts));
}

bool Plugin::PhysicsBubbleIsEmpty() const {
  VLOG(1) << __FUNCTION__;
  VLOG_AND_RETURN(1, bubble_->empty());
}

Displacement<World> Plugin::BubbleDisplacementCorrection(
    Position<World> const& sun_world_position) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(sun_world_position);
  VLOG_AND_RETURN(1, bubble_->DisplacementCorrection(BarycentricToWorldSun(),
                                                     *sun_,
                                                     sun_world_position));
}

Velocity<World> Plugin::BubbleVelocityCorrection(
    Index const reference_body_index) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(reference_body_index);
  Celestial const& reference_body =
      *FindOrDie(celestials_, reference_body_index);
  VLOG_AND_RETURN(1, bubble_->VelocityCorrection(BarycentricToWorldSun(),
                                                 reference_body));
}

FrameField<World> Plugin::Navball(
    Position<World> const& sun_world_position) const {
  auto const to_world =
      OrthogonalMap<WorldSun, World>::Identity() * BarycentricToWorldSun();
  ephemeris_->Prolong(current_time_);
  auto const positions_from_world =
      AffineMap<World, Barycentric, Length, OrthogonalMap>(
          sun_world_position,
          sun_->current_position(current_time_),
          to_world.Inverse());
  return [this, positions_from_world, to_world](
      Position<World> const& q) -> Rotation<World, World> {
    // KSP's navball has x west, y up, z south.
    // we want x north, y west, z up.
    auto const orthogonal_map = to_world *
        plotting_frame_->FromThisFrameAtTime(current_time_).orthogonal_map() *
        Permutation<World, Navigation>(
            Permutation<World, Navigation>::XZY).Forget() *
        Rotation<World, World>(π / 2 * Radian,
                               Bivector<double, World>({0, 1, 0})).Forget();
    CHECK(orthogonal_map.Determinant().Positive());
    return orthogonal_map.rotation();
  };
}

Vector<double, World> Plugin::VesselTangent(GUID const& vessel_guid) const {
  return FromVesselFrenetFrame(*find_vessel_by_guid_or_die(vessel_guid),
                               Vector<double, Frenet<Navigation>>({1, 0, 0}));
}

Vector<double, World> Plugin::VesselNormal(GUID const& vessel_guid) const {
  return FromVesselFrenetFrame(*find_vessel_by_guid_or_die(vessel_guid),
                               Vector<double, Frenet<Navigation>>({0, 1, 0}));
}

Vector<double, World> Plugin::VesselBinormal(GUID const& vessel_guid) const {
  return FromVesselFrenetFrame(*find_vessel_by_guid_or_die(vessel_guid),
                               Vector<double, Frenet<Navigation>>({0, 0, 1}));
}

OrthogonalMap<Barycentric, WorldSun> Plugin::BarycentricToWorldSun() const {
  return kSunLookingGlass.Inverse().Forget() * PlanetariumRotation().Forget();
}

Instant Plugin::CurrentTime() const {
  return current_time_;
}

void Plugin::WriteToMessage(
    not_null<serialization::Plugin*> const message) const {
  LOG(INFO) << __FUNCTION__;
  CHECK(!initializing_);
  ephemeris_->Prolong(current_time_);
  std::map<not_null<Celestial const*>, Index const> celestial_to_index;
  for (auto const& index_celestial : celestials_) {
    celestial_to_index.emplace(index_celestial.second.get(),
                               index_celestial.first);
  }
  for (auto const& index_celestial : celestials_) {
    Index const index = index_celestial.first;
    not_null<Celestial const*> const celestial = index_celestial.second.get();
    auto* const celestial_message = message->add_celestial();
    celestial_message->set_index(index);
    if (celestial->has_parent()) {
      Index const parent_index =
          FindOrDie(celestial_to_index, celestial->parent());
      celestial_message->set_parent_index(parent_index);
    }
  }
  std::map<not_null<Vessel const*>, GUID const> vessel_to_guid;
  for (auto const& guid_vessel : vessels_) {
    std::string const& guid = guid_vessel.first;
    not_null<Vessel*> const vessel = guid_vessel.second.get();
    vessel_to_guid.emplace(vessel, guid);
    auto* const vessel_message = message->add_vessel();
    vessel_message->set_guid(guid);
    vessel->WriteToMessage(vessel_message->mutable_vessel());
    Index const parent_index = FindOrDie(celestial_to_index, vessel->parent());
    vessel_message->set_parent_index(parent_index);
    vessel_message->set_dirty(is_dirty(vessel));
  }

  ephemeris_->WriteToMessage(message->mutable_ephemeris());
  prolongation_integrator_.WriteToMessage(
      message->mutable_prolongation_integrator());
  prediction_integrator_.WriteToMessage(
      message->mutable_prediction_integrator());

  bubble_->WriteToMessage(
      [&vessel_to_guid](not_null<Vessel const*> const vessel) -> GUID {
        return FindOrDie(vessel_to_guid, vessel);
      },
      message->mutable_bubble());

  planetarium_rotation_.WriteToMessage(message->mutable_planetarium_rotation());
  current_time_.WriteToMessage(message->mutable_current_time());
  Index const sun_index = FindOrDie(celestial_to_index, sun_);
  message->set_sun_index(sun_index);
  plotting_frame_->WriteToMessage(message->mutable_plotting_frame());
  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

not_null<std::unique_ptr<Plugin>> Plugin::ReadFromMessage(
    serialization::Plugin const& message) {
  LOG(INFO) << __FUNCTION__;
  bool const is_pre_bourbaki = message.pre_bourbaki_celestial_size() > 0;
  std::unique_ptr<Ephemeris<Barycentric>> ephemeris;
  IndexToOwnedCelestial celestials;

  if (is_pre_bourbaki) {
    ephemeris = Ephemeris<Barycentric>::ReadFromPreBourbakiMessages(
        message.pre_bourbaki_celestial(),
        McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
        kStep,
        kFittingTolerance);
    ReadCelestialsFromMessages(*ephemeris,
                               message.pre_bourbaki_celestial(),
                               &celestials);
  } else {
    ephemeris = Ephemeris<Barycentric>::ReadFromMessage(message.ephemeris());
    ReadCelestialsFromMessages(*ephemeris,
                               message.celestial(),
                               &celestials);
  }

  GUIDToOwnedVessel vessels;
  std::set<not_null<Vessel*>> dirty_vessels;
  for (auto const& vessel_message : message.vessel()) {
    not_null<Celestial const*> const parent =
        FindOrDie(celestials, vessel_message.parent_index()).get();
    not_null<std::unique_ptr<Vessel>> vessel = Vessel::ReadFromMessage(
                                                   vessel_message.vessel(),
                                                   ephemeris.get(),
                                                   parent);
    if (vessel_message.dirty()) {
      dirty_vessels.emplace(vessel.get());
    }
    auto const inserted =
        vessels.emplace(vessel_message.guid(), std::move(vessel));
    CHECK(inserted.second);
  }
  not_null<std::unique_ptr<PhysicsBubble>> bubble =
      PhysicsBubble::ReadFromMessage(
          [&vessels](GUID guid) -> not_null<Vessel*> {
            return FindOrDie(vessels, guid).get();
          },
          message.bubble());
  auto const& prolongation_integrator =
      is_pre_bourbaki ?
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>() :
          AdaptiveStepSizeIntegrator<NewtonianMotionEquation>::ReadFromMessage(
              message.prolongation_integrator());
  auto const& prediction_integrator =
      is_pre_bourbaki ?
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>() :
          AdaptiveStepSizeIntegrator<NewtonianMotionEquation>::ReadFromMessage(
              message.prediction_integrator());

  Instant const current_time = Instant::ReadFromMessage(message.current_time());
  Instant history_time = current_time;
  for (auto const& pair : vessels) {
    auto const& vessel = pair.second;
    if (vessel->is_synchronized()) {
      history_time = vessel->history().last().time();
      break;
    }
  }

  // Can't use |make_unique| here without implementation-dependent friendships.
  auto plugin = std::unique_ptr<Plugin>(
      new Plugin(std::move(vessels),
                 std::move(celestials),
                 std::move(dirty_vessels),
                 std::move(bubble),
                 std::move(ephemeris),
                 prolongation_integrator,
                 prediction_integrator,
                 Angle::ReadFromMessage(message.planetarium_rotation()),
                 current_time,
                 history_time,
                 message.sun_index()));
  std::unique_ptr<NavigationFrame> plotting_frame =
      NavigationFrame::ReadFromMessage(plugin->ephemeris_.get(),
                                       message.plotting_frame());
  if (plotting_frame == nullptr) {
    // In the pre-Brouwer compatibility case you get a plotting frame centered
    // on the Sun.
    plugin->SetPlottingFrame(
        plugin->NewBodyCentredNonRotatingNavigationFrame(message.sun_index()));
  } else {
    plugin->SetPlottingFrame(std::move(plotting_frame));
  }
  return std::move(plugin);
}

Plugin::Plugin(GUIDToOwnedVessel vessels,
               IndexToOwnedCelestial celestials,
               std::set<not_null<Vessel*>> dirty_vessels,
               not_null<std::unique_ptr<PhysicsBubble>> bubble,
               std::unique_ptr<Ephemeris<Barycentric>> ephemeris,
               AdaptiveStepSizeIntegrator<
                   NewtonianMotionEquation> const& prolongation_integrator,
               AdaptiveStepSizeIntegrator<
                   NewtonianMotionEquation> const& prediction_integrator,
               Angle planetarium_rotation,
               Instant current_time,
               Instant history_time,
               Index sun_index)
    : vessels_(std::move(vessels)),
      celestials_(std::move(celestials)),
      dirty_vessels_(std::move(dirty_vessels)),
      bubble_(std::move(bubble)),
      ephemeris_(std::move(ephemeris)),
      history_integrator_(ephemeris_->planetary_integrator()),
      prolongation_integrator_(prolongation_integrator),
      prediction_integrator_(prediction_integrator),
      planetarium_rotation_(planetarium_rotation),
      current_time_(current_time),
      history_time_(history_time),
      sun_(FindOrDie(celestials_, sun_index).get()) {
  for (auto const& guid_vessel : vessels_) {
    auto const& vessel = guid_vessel.second;
    kept_vessels_.emplace(vessel.get());
    if (!vessel->is_synchronized()) {
      unsynchronized_vessels_.emplace(vessel.get());
    }
  }
  initializing_.Flop();
}


void Plugin::InitializeEphemerisAndSetCelestialTrajectories() {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<Barycentric>> initial_state;
  for (auto& pair : absolute_initialization_->bodies) {
    auto& body = pair.second;
    bodies.emplace_back(std::move(body));
  }
  for (auto const& state : absolute_initialization_->initial_state) {
    initial_state.emplace_back(state.second);
  }
  absolute_initialization_ = std::experimental::nullopt;
  ephemeris_ = std::make_unique<Ephemeris<Barycentric>>(
      std::move(bodies),
      initial_state,
      current_time_,
      history_integrator_,
      kStep,
      kFittingTolerance);
  for (auto const& pair : celestials_) {
    auto& celestial = *pair.second;
    celestial.set_trajectory(ephemeris_->trajectory(celestial.body()));
  }

  // This would use NewBodyCentredNonRotatingNavigationFrame, but we don't have
  // the sun's index at hand.
  // TODO(egg): maybe these functions should take |Celestial*|s, and we should
  // then export |FindOrDie(celestials_, _)|.
  SetPlottingFrame(
      make_not_null_unique<
          BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>>(
          ephemeris_.get(),
          sun_->body()));
}

not_null<std::unique_ptr<Vessel>> const& Plugin::find_vessel_by_guid_or_die(
    GUID const& vessel_guid) const {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(vessel_guid);
  VLOG_AND_RETURN(1, FindOrDie(vessels_, vessel_guid));
}

bool Plugin::has_dirty_vessels() const {
  return !dirty_vessels_.empty();
}

bool Plugin::has_unsynchronized_vessels() const {
  return !unsynchronized_vessels_.empty();
}

bool Plugin::is_dirty(not_null<Vessel*> const vessel) const {
  return dirty_vessels_.count(vessel) > 0;
}

// The map between the vector spaces of |Barycentric| and |AliceSun| at
// |current_time_|.
Rotation<Barycentric, AliceSun> Plugin::PlanetariumRotation() const {
  return Rotation<Barycentric, AliceSun>(
      planetarium_rotation_,
      Bivector<double, Barycentric>({0, 0, -1}));
}

void Plugin::CleanUpVessels() {
  VLOG(1) <<  __FUNCTION__;
  // Remove the vessels which were not updated since last time.
  for (auto it = vessels_.cbegin(); it != vessels_.cend();) {
    // While we're going over the vessels, check invariants.
    CheckVesselInvariants(it);
    not_null<Vessel*> const vessel = it->second.get();
    // Now do the cleanup.
    if (kept_vessels_.erase(vessel)) {
      ++it;
    } else {
      LOG(INFO) << "Removing vessel with GUID " << it->first;
      // Since we are going to delete the vessel, we must remove it from
      // |new_vessels| if it's there.
      if (unsynchronized_vessels_.erase(vessel)) {
        LOG(INFO) << "Vessel had not been synchronized";
      }
      if (dirty_vessels_.erase(vessel)) {
        LOG(INFO) << "Vessel was dirty";
      }
      // |std::map::erase| invalidates its parameter so we post-increment.
      vessels_.erase(it++);
    }
  }
}

void Plugin::CheckVesselInvariants(
    GUIDToOwnedVessel::const_iterator const it) const {
  not_null<std::unique_ptr<Vessel>> const& vessel = it->second;
  CHECK(vessel->is_initialized()) << "Vessel with GUID " << it->first
                                  << " was not given an initial state";
  // TODO(egg): At the moment, if a vessel is inserted when
  // |current_time_ == HistoryTime()| (that only happens before the first call
  // to |AdvanceTime|) its first step is unsynchronized. This is convenient to
  // test code paths, but it means the invariant is GE, rather than GT.
  CHECK_GE(vessel->prolongation().last().time(), history_time_);
  if (unsynchronized_vessels_.count(vessel.get()) > 0) {
    CHECK(!vessel->is_synchronized());
  } else {
    CHECK(vessel->is_synchronized());
    CHECK_EQ(vessel->history().last().time(), history_time_);
  }
}

Plugin::Trajectories Plugin::SynchronizedHistories() const {
  Trajectories trajectories;
  // NOTE(egg): This may be too large, vessels that are not new and in the
  // physics bubble or dirty will not be added.
  trajectories.reserve(vessels_.size() - unsynchronized_vessels_.size());
  for (auto const& pair : vessels_) {
    not_null<Vessel*> const vessel = pair.second.get();
    if (vessel->is_synchronized() &&
        !bubble_->contains(vessel) &&
        !is_dirty(vessel)) {
      trajectories.push_back(vessel->mutable_history());
    }
  }
  return trajectories;
}

void Plugin::EvolveHistories(Instant const& t, Trajectories const& histories) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(t);
  // Integration with a constant step.{
  VLOG(1) << "Starting the evolution of the histories" << '\n'
          << "from : " << history_time_;
  ephemeris_->FlowWithFixedStep(
      histories,
      Ephemeris<Barycentric>::kNoIntrinsicAccelerations,
      Δt_,
      t);
  history_time_ = histories.front()->last().time();
  CHECK_GE(history_time_, current_time_);
  VLOG(1) << "Evolved the histories" << '\n'
          << "to   : " << history_time_;
}

void Plugin::SynchronizeNewVesselsAndCleanDirtyVessels() {
  VLOG(1) << __FUNCTION__;
  Trajectories trajectories;
  Ephemeris<Barycentric>::IntrinsicAccelerations intrinsic_accelerations;

  trajectories.reserve(unsynchronized_vessels_.size() + bubble_->count());
  intrinsic_accelerations.reserve(
      unsynchronized_vessels_.size() + bubble_->count());
  for (not_null<Vessel*> const vessel : unsynchronized_vessels_) {
    if (!bubble_->contains(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
      intrinsic_accelerations.push_back(
          Ephemeris<Barycentric>::kNoIntrinsicAcceleration);
    }
  }
  for (not_null<Vessel*> const vessel : dirty_vessels_) {
    if (!bubble_->contains(vessel) && vessel->is_synchronized()) {
      trajectories.push_back(vessel->mutable_prolongation());
      intrinsic_accelerations.push_back(
          Ephemeris<Barycentric>::kNoIntrinsicAcceleration);
    }
  }
  if (!bubble_->empty()) {
    trajectories.push_back(bubble_->mutable_centre_of_mass_trajectory());
    intrinsic_accelerations.push_back(
        bubble_->centre_of_mass_intrinsic_acceleration());
  }
  VLOG(1) << "Starting the synchronization of the new vessels"
          << (bubble_->empty() ? "" : " and of the bubble");
  for (int i = 0; i < trajectories.size(); ++i) {
    auto const& trajectory = trajectories[i];
    auto const& intrinsic_acceleration = intrinsic_accelerations[i];
    ephemeris_->FlowWithAdaptiveStep(trajectory,
                                     intrinsic_acceleration,
                                     prolongation_length_tolerance_,
                                     prolongation_speed_tolerance_,
                                     prolongation_integrator_,
                                     history_time_);
  }
  if (!bubble_->empty()) {
    SynchronizeBubbleHistories();
  }
  for (not_null<Vessel*> const vessel : unsynchronized_vessels_) {
    CHECK(!bubble_->contains(vessel));
    vessel->CreateHistoryAndForkProlongation(
        history_time_,
        vessel->prolongation().last().degrees_of_freedom());
    dirty_vessels_.erase(vessel);
  }
  unsynchronized_vessels_.clear();
  for (not_null<Vessel*> const vessel : dirty_vessels_) {
    CHECK(!bubble_->contains(vessel));
    vessel->mutable_history()->Append(
        history_time_,
        vessel->prolongation().last().degrees_of_freedom());
  }
  dirty_vessels_.clear();
  VLOG(1) << "Synchronized the new vessels"
          << (bubble_->empty() ? "" : " and the bubble");
}

void Plugin::SynchronizeBubbleHistories() {
  VLOG(1) << __FUNCTION__;
  DegreesOfFreedom<Barycentric> const& centre_of_mass =
      bubble_->centre_of_mass_trajectory().last().degrees_of_freedom();
  for (not_null<Vessel*> const vessel : bubble_->vessels()) {
    RelativeDegreesOfFreedom<Barycentric> const& from_centre_of_mass =
        bubble_->from_centre_of_mass(vessel);
    if (vessel->is_synchronized()) {
      vessel->mutable_history()->Append(
          history_time_,
          centre_of_mass + from_centre_of_mass);
    } else {
      vessel->CreateHistoryAndForkProlongation(
          history_time_,
          centre_of_mass + from_centre_of_mass);
      CHECK(unsynchronized_vessels_.erase(vessel));
    }
    CHECK(dirty_vessels_.erase(vessel));
  }
}

void Plugin::ResetProlongations() {
  VLOG(1) << __FUNCTION__;
  for (auto const& pair : vessels_) {
    not_null<std::unique_ptr<Vessel>> const& vessel = pair.second;
    vessel->ResetProlongation(history_time_);
  }
  VLOG(1) << "Prolongations have been reset";
}

void Plugin::EvolveProlongationsAndBubble(Instant const& t) {
  VLOG(1) << __FUNCTION__ << '\n' << NAMED(t);
  Trajectories trajectories;
  Ephemeris<Barycentric>::IntrinsicAccelerations intrinsic_accelerations;

  std::size_t const number_of_vessels_not_in_physics_bubble =
      vessels_.size() - bubble_->number_of_vessels();
  trajectories.reserve(number_of_vessels_not_in_physics_bubble +
                       bubble_->count());
  intrinsic_accelerations.reserve(number_of_vessels_not_in_physics_bubble +
                                  bubble_->count());
  for (auto const& pair : vessels_) {
    not_null<Vessel*> const vessel = pair.second.get();
    if (!bubble_->contains(vessel)) {
      trajectories.push_back(vessel->mutable_prolongation());
      intrinsic_accelerations.push_back(
          Ephemeris<Barycentric>::kNoIntrinsicAcceleration);
    }
  }
  if (!bubble_->empty()) {
    trajectories.push_back(bubble_->mutable_centre_of_mass_trajectory());
    intrinsic_accelerations.push_back(
        bubble_->centre_of_mass_intrinsic_acceleration());
  }
  VLOG(1) << "Evolving prolongations"
          << (bubble_->empty() ? "" : " and bubble") << '\n'
          << "from : " << trajectories.front()->last().time() << '\n'
          << "to   : " << t;
  for (int i = 0; i < trajectories.size(); ++i) {
    auto const& trajectory = trajectories[i];
    auto const& intrinsic_acceleration = intrinsic_accelerations[i];
    ephemeris_->FlowWithAdaptiveStep(trajectory,
                                     intrinsic_acceleration,
                                     prolongation_length_tolerance_,
                                     prolongation_speed_tolerance_,
                                     prolongation_integrator_,
                                     t);
  }
  if (!bubble_->empty()) {
    DegreesOfFreedom<Barycentric> const& centre_of_mass =
        bubble_->centre_of_mass_trajectory().last().degrees_of_freedom();
    for (not_null<Vessel*> vessel : bubble_->vessels()) {
      RelativeDegreesOfFreedom<Barycentric> const& from_centre_of_mass =
          bubble_->from_centre_of_mass(vessel);
      vessel->mutable_prolongation()->Append(
          t,
          centre_of_mass + from_centre_of_mass);
    }
  }
}

Vector<double, World> Plugin::FromVesselFrenetFrame(
    Vessel const& vessel,
    Vector<double, Frenet<Navigation>> const& vector) const {
  auto const& last = vessel.prolongation().last();
  Instant const& time = last.time();
  DegreesOfFreedom<Barycentric> const& degrees_of_freedom =
      last.degrees_of_freedom();
  auto const from_frenet_frame_to_navigation_frame =
      plotting_frame_->FrenetFrame(
          time,
          plotting_frame_->ToThisFrameAtTime(time)(degrees_of_freedom));

  // The given |vector| in the Frenet frame of the vessel's free-falling
  // trajectory in the given |navigation_frame|, converted to |WorldSun|
  // coordinates.
  return Identity<WorldSun, World>()(
      BarycentricToWorldSun()(
          plotting_frame_->FromThisFrameAtTime(time).orthogonal_map()(
              from_frenet_frame_to_navigation_frame(vector))));
}

template<typename T>
void Plugin::ReadCelestialsFromMessages(
  Ephemeris<Barycentric> const& ephemeris,
  google::protobuf::RepeatedPtrField<T> const& celestial_messages,
  not_null<IndexToOwnedCelestial*> const celestials) {
  auto const& bodies = ephemeris.bodies();
  auto bodies_it = bodies.begin();
  for (auto const& celestial_message : celestial_messages) {
    auto const inserted =
      celestials->emplace(celestial_message.index(),
                          make_not_null_unique<Celestial>(*bodies_it));
    CHECK(inserted.second);
    inserted.first->second->set_trajectory(ephemeris.trajectory(*bodies_it));
    ++bodies_it;
  }
  CHECK_EQ(bodies.end() - bodies.begin(), bodies_it - bodies.begin());
  for (auto const& celestial_message : celestial_messages) {
    if (celestial_message.has_parent_index()) {
      not_null<std::unique_ptr<Celestial>> const& celestial =
        FindOrDie(*celestials, celestial_message.index());
      not_null<Celestial const*> const parent =
        FindOrDie(*celestials, celestial_message.parent_index()).get();
      celestial->set_parent(parent);
    }
  }
}

}  // namespace ksp_plugin
}  // namespace principia
