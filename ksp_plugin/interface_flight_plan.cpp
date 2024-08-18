#include "ksp_plugin/interface.hpp"

#include <vector>
#include <utility>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "glog/logging.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/flight_plan_optimization_driver.hpp"
#include "ksp_plugin/flight_plan_optimizer.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/renderer.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/barycentric_rotating_reference_frame.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_flight_plan_optimization_driver;
using namespace principia::ksp_plugin::_flight_plan_optimizer;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_iterators;
using namespace principia::ksp_plugin::_renderer;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_barycentric_rotating_reference_frame;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

namespace {

NavigationManœuvre::Burn FromInterfaceBurn(Plugin const& plugin,
                                           Burn const& burn) {
  NavigationManœuvre::Intensity intensity;
  intensity.Δv = FromXYZ<Velocity<Frenet<NavigationFrame>>>(burn.delta_v);
  NavigationManœuvre::Timing timing;
  timing.initial_time = FromGameTime(plugin, burn.initial_time);
  return {intensity,
          timing,
          burn.thrust_in_kilonewtons * Kilo(Newton),
          burn.specific_impulse_in_seconds_g0 * Second * StandardGravity,
          NewNavigationFrame(plugin, burn.frame),
          burn.is_inertially_fixed};
}

Burn GetBurn(Plugin const& plugin,
             NavigationManœuvre const& manœuvre) {
  // When building the parameters, make sure that the "optional" fields get a
  // deterministic default.
  NavigationFrameParameters parameters;
  parameters.centre_index = -1;
  parameters.primary_index = -1;
  parameters.secondary_index = -1;

  int number_of_subclasses = 0;

  {
    auto const* barycentric_rotating_reference_frame = dynamic_cast<
        BarycentricRotatingReferenceFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (barycentric_rotating_reference_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension = serialization::BarycentricRotatingReferenceFrame::
          kExtensionFieldNumber;
      // A barycentric rotating reference frame can have multiple secondaries
      // since ابن الهيثم and multiple primaries since 伊藤, but it has not been
      // possible to set a manœuvre frame to barycentric since Haar, so this
      // cannot happen in this compatibility path.
      CHECK_EQ(barycentric_rotating_reference_frame->primaries().size(), 1);
      parameters.primary_index = plugin.CelestialIndexOfBody(
          *barycentric_rotating_reference_frame->primaries().front());
      CHECK_EQ(barycentric_rotating_reference_frame->secondaries().size(), 1);
      parameters.secondary_index = plugin.CelestialIndexOfBody(
          *barycentric_rotating_reference_frame->secondaries().front());
    }
  }

  {
    auto const* body_centred_body_direction_reference_frame = dynamic_cast<
        BodyCentredBodyDirectionReferenceFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (body_centred_body_direction_reference_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension = serialization::
          BodyCentredBodyDirectionReferenceFrame::kExtensionFieldNumber;
      parameters.primary_index = plugin.CelestialIndexOfBody(
          *body_centred_body_direction_reference_frame->primary());
      parameters.secondary_index = plugin.CelestialIndexOfBody(
          *body_centred_body_direction_reference_frame->secondary());
    }
  }

  {
    auto const* body_centred_non_rotating_reference_frame = dynamic_cast<
        BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (body_centred_non_rotating_reference_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension = serialization::
          BodyCentredNonRotatingReferenceFrame::kExtensionFieldNumber;
      parameters.centre_index = plugin.CelestialIndexOfBody(
          *body_centred_non_rotating_reference_frame->centre());
    }
  }

  {
    auto const* body_surface_reference_frame = dynamic_cast<
        BodySurfaceReferenceFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (body_surface_reference_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension =
          serialization::BodySurfaceReferenceFrame::kExtensionFieldNumber;
      parameters.centre_index =
          plugin.CelestialIndexOfBody(*body_surface_reference_frame->centre());
    }
  }

  CHECK_EQ(number_of_subclasses, 1) << "Could not construct frame parameters";

  return {manœuvre.thrust() / Kilo(Newton),
          manœuvre.specific_impulse() / (Second * StandardGravity),
          parameters,
          ToGameTime(plugin, manœuvre.initial_time()),
          ToXYZ(manœuvre.Δv()),
          manœuvre.is_inertially_fixed()};
}

NavigationManoeuvre* ToNewInterfaceNavigationManoeuvre(
    Plugin const& plugin,
    NavigationManœuvre const& manœuvre) {
  return new NavigationManoeuvre{
      .burn = GetBurn(plugin, manœuvre),
      .initial_mass_in_tonnes = manœuvre.initial_mass() / Tonne,
      .final_mass_in_tonnes = manœuvre.final_mass() / Tonne,
      .mass_flow = manœuvre.mass_flow() / (Kilogram / Second),
      .duration = manœuvre.duration() / Second,
      .final_time = ToGameTime(plugin, manœuvre.final_time()),
      .time_of_half_delta_v = ToGameTime(plugin, manœuvre.time_of_half_Δv()),
      .time_to_half_delta_v = manœuvre.time_to_half_Δv() / Second};
}

}  // namespace

int __cdecl principia__FlightPlanCount(Plugin const* const plugin,
                                       char const* const vessel_guid) {
  journal::Method<journal::FlightPlanCount> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(plugin->GetVessel(vessel_guid)->flight_plan_count());
}

void __cdecl principia__FlightPlanCreate(Plugin const* const plugin,
                                         char const* const vessel_guid,
                                         double const final_time,
                                         double const mass_in_tonnes) {
  journal::Method<journal::FlightPlanCreate> m({plugin,
                                                vessel_guid,
                                                final_time,
                                                mass_in_tonnes});
  CHECK_NOTNULL(plugin);
  plugin->CreateFlightPlan(vessel_guid,
                           FromGameTime(*plugin, final_time),
                           mass_in_tonnes * Tonne);
  return m.Return();
}

void __cdecl principia__FlightPlanDelete(Plugin const* const plugin,
                                         char const* const vessel_guid) {
  journal::Method<journal::FlightPlanDelete> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  plugin->GetVessel(vessel_guid)->DeleteFlightPlan();
  return m.Return();
}

void __cdecl principia__FlightPlanDuplicate(Plugin const* const plugin,
                                            char const* const vessel_guid) {
  journal::Method<journal::FlightPlanDuplicate> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  plugin->GetVessel(vessel_guid)->DuplicateFlightPlan();
  return m.Return();
}

bool __cdecl principia__FlightPlanExists(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanExists> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(plugin->GetVessel(vessel_guid)->has_flight_plan());
}

double __cdecl principia__FlightPlanGetActualFinalTime(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetActualFinalTime> m(
      {plugin, vessel_guid});
  return m.Return(
      ToGameTime(*plugin,
                 GetFlightPlan(*plugin, vessel_guid).actual_final_time()));
}

FlightPlanAdaptiveStepParameters __cdecl
principia__FlightPlanGetAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetAdaptiveStepParameters> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  auto const& flight_plan = GetFlightPlan(*plugin, vessel_guid);
  return m.Return(ToFlightPlanAdaptiveStepParameters(
                      flight_plan.adaptive_step_parameters(),
                      flight_plan.generalized_adaptive_step_parameters()));
}

Status* __cdecl principia__FlightPlanGetAnomalousStatus(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetAnomalousStatus> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(
      ToNewStatus(GetFlightPlan(*plugin, vessel_guid).anomalous_status()));
}

OrbitAnalysis* __cdecl principia__FlightPlanGetCoastAnalysis(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const* const revolutions_per_cycle,
    int const* const days_per_cycle,
    int const ground_track_revolution,
    int const index) {
  journal::Method<journal::FlightPlanGetCoastAnalysis> m(
      {plugin,
       vessel_guid,
       revolutions_per_cycle,
       days_per_cycle,
       ground_track_revolution,
       index});
  CHECK_NOTNULL(plugin);
  auto& flight_plan = GetFlightPlan(*plugin, vessel_guid);
  auto const analysis =
      NewOrbitAnalysis(flight_plan.analysis(index),
                       *plugin,
                       /*revolutions_per_cycle=*/revolutions_per_cycle,
                       /*days_per_cycle=*/days_per_cycle,
                       /*ground_track_revolution=*/ground_track_revolution);
  analysis->progress_of_next_analysis = flight_plan.progress_of_analysis(index);
  return m.Return(analysis);
}

double __cdecl principia__FlightPlanGetDesiredFinalTime(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetDesiredFinalTime> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(
      ToGameTime(*plugin,
                 GetFlightPlan(*plugin, vessel_guid).desired_final_time()));
}

XYZ __cdecl principia__FlightPlanGetGuidance(Plugin const* const plugin,
                                             char const* const vessel_guid,
                                             int const index) {
  journal::Method<journal::FlightPlanGetGuidance> m(
      {plugin, vessel_guid, index});
  CHECK_NOTNULL(plugin);
  auto const& manœuvre = GetFlightPlan(*plugin, vessel_guid).GetManœuvre(index);
  Vector<double, World> result;
  if (manœuvre.is_inertially_fixed()) {
    result = plugin->renderer().BarycentricToWorld(
                 plugin->PlanetariumRotation())(manœuvre.InertialDirection());
  } else {
    result = (plugin->CurrentTime() < manœuvre.initial_time()
                  ? plugin->renderer().BarycentricToWorld(
                        plugin->PlanetariumRotation()) *
                        manœuvre.FrenetFrame()
                  : plugin->renderer().FrenetToWorld(
                        *plugin->GetVessel(vessel_guid),
                        *manœuvre.frame(),
                        plugin->PlanetariumRotation()))(manœuvre.direction());
  }
  return m.Return(ToXYZ(result));
}

double __cdecl principia__FlightPlanGetInitialTime(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetInitialTime> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(
      ToGameTime(*plugin,
                 GetFlightPlan(*plugin, vessel_guid).initial_time()));
}

NavigationManoeuvre* __cdecl principia__FlightPlanGetManoeuvre(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::FlightPlanGetManoeuvre> m({plugin,
                                                      vessel_guid,
                                                      index});
  CHECK_NOTNULL(plugin);
  return m.Return(ToNewInterfaceNavigationManoeuvre(
                      *plugin,
                      GetFlightPlan(*plugin, vessel_guid).GetManœuvre(index)));
}

NavigationManoeuvreFrenetTrihedron __cdecl
principia__FlightPlanGetManoeuvreFrenetTrihedron(Plugin const* const plugin,
                                                 char const* const vessel_guid,
                                                 int const index) {
  journal::Method<journal::FlightPlanGetManoeuvreFrenetTrihedron> m(
      {plugin, vessel_guid, index});
  CHECK_NOTNULL(plugin);

  NavigationManœuvre const& manœuvre =
      GetFlightPlan(*plugin, vessel_guid).GetManœuvre(index);
  OrthogonalMap<Frenet<Navigation>, World> const frenet_to_plotted_world =
      plugin->renderer().FrenetToWorld(plugin->CurrentTime(),
                                       manœuvre,
                                       plugin->PlanetariumRotation());
  NavigationManoeuvreFrenetTrihedron result;
  result.tangent = ToXYZ(
      frenet_to_plotted_world(Vector<double, Frenet<Navigation>>({1, 0, 0})));
  result.normal = ToXYZ(
      frenet_to_plotted_world(Vector<double, Frenet<Navigation>>({0, 1, 0})));
  result.binormal = ToXYZ(
      frenet_to_plotted_world(Vector<double, Frenet<Navigation>>({0, 0, 1})));

  return m.Return(result);
}

XYZ __cdecl principia__FlightPlanGetManoeuvreInitialPlottedVelocity(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::FlightPlanGetManoeuvreInitialPlottedVelocity> m(
      {plugin, vessel_guid, index});
  CHECK_NOTNULL(plugin);

  auto const& [t, dof] =
      GetFlightPlan(*plugin, vessel_guid).GetSegment(2 * index)->back();
  Velocity<Navigation> const v =
      plugin->renderer().BarycentricToPlotting(t)(dof).velocity();
  return m.Return(ToXYZ(plugin->renderer().PlottingToWorld(
      plugin->CurrentTime(), plugin->PlanetariumRotation())(v)));
}

Status* __cdecl principia__FlightPlanInsert(Plugin const* const plugin,
                                            char const* const vessel_guid,
                                            Burn const& burn,
                                            int const index) {
  journal::Method<journal::FlightPlanInsert> m(
      {plugin, vessel_guid, burn, index});
  CHECK_NOTNULL(plugin);
  auto const status = GetFlightPlan(*plugin, vessel_guid)
                          .Insert(FromInterfaceBurn(*plugin, burn), index);
  plugin->ExtendPredictionForFlightPlan(vessel_guid);
  return m.Return(ToNewStatus(status));
}

int __cdecl principia__FlightPlanNumberOfAnomalousManoeuvres(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfAnomalousManoeuvres> m(
      {plugin,
       vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(GetFlightPlan(*plugin, vessel_guid).
                      number_of_anomalous_manœuvres());
}

int __cdecl principia__FlightPlanNumberOfManoeuvres(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfManoeuvres> m({plugin,
                                                            vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(GetFlightPlan(*plugin, vessel_guid).number_of_manœuvres());
}

int __cdecl principia__FlightPlanNumberOfSegments(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfSegments> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(GetFlightPlan(*plugin, vessel_guid).number_of_segments());
}

int __cdecl principia__FlightPlanOptimizationDriverInProgress(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanOptimizationDriverInProgress> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  auto& vessel = *plugin->GetVessel(vessel_guid);
  auto const maybe_parameters = vessel.FlightPlanOptimizationDriverInProgress();
  if (maybe_parameters.has_value()) {
    return m.Return(maybe_parameters->index);
  } else {
    return m.Return(-1);
  }
}

void __cdecl principia__FlightPlanOptimizationDriverMake(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const distance,
    double const* const inclination_in_degrees,
    int const celestial_index,
    NavigationFrameParameters const& navigation_frame_parameters) {
  journal::Method<journal::FlightPlanOptimizationDriverMake> m(
      {plugin,
       vessel_guid,
       distance,
       inclination_in_degrees,
       celestial_index,
       navigation_frame_parameters});
  CHECK_NOTNULL(plugin);
  auto& vessel = *plugin->GetVessel(vessel_guid);
  const auto& celestial = plugin->GetCelestial(celestial_index);

  std::vector<FlightPlanOptimizer::MetricFactory> factories = {
      FlightPlanOptimizer::ForCelestialDistance(
          /*celestial=*/&celestial,
          /*target_distance=*/distance * Metre),
      FlightPlanOptimizer::ForΔv()};
  std::vector<double> weights = {1, 1e3};
  if (inclination_in_degrees != nullptr) {
    factories.push_back(FlightPlanOptimizer::ForInclination(
        &celestial,
        [plugin, navigation_frame_parameters]() {
          return NewNavigationFrame(*plugin, navigation_frame_parameters);
        },
        *inclination_in_degrees * Degree));
    weights.push_back(1e6);
  }

  vessel.MakeFlightPlanOptimizationDriver(
      FlightPlanOptimizer::LinearCombination(factories, weights));

  return m.Return();
}

void __cdecl principia__FlightPlanOptimizationDriverStart(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const manœuvre_index) {
  journal::Method<journal::FlightPlanOptimizationDriverStart> m(
      {plugin,
       vessel_guid,
       manœuvre_index});
  CHECK_NOTNULL(plugin);
  auto& vessel = *plugin->GetVessel(vessel_guid);

  vessel.StartFlightPlanOptimizationDriver(
      {.index = manœuvre_index,
       .Δv_tolerance = 1 * Micro(Metre) / Second});

  return m.Return();
}

Status* __cdecl principia__FlightPlanRebase(Plugin const* const plugin,
                                         char const* const vessel_guid,
                                         double const mass_in_tonnes) {
  journal::Method<journal::FlightPlanRebase> m(
      {plugin, vessel_guid, mass_in_tonnes});
  CHECK_NOTNULL(plugin);
  auto const status =
      plugin->GetVessel(vessel_guid)->RebaseFlightPlan(mass_in_tonnes * Tonne);
  plugin->ExtendPredictionForFlightPlan(vessel_guid);
  return m.Return(ToNewStatus(status));
}

Status* __cdecl principia__FlightPlanRemove(Plugin const* const plugin,
                                            char const* const vessel_guid,
                                            int const index) {
  journal::Method<journal::FlightPlanRemove> m({plugin, vessel_guid, index});
  CHECK_NOTNULL(plugin);
  auto const status = GetFlightPlan(*plugin, vessel_guid).Remove(index);
  plugin->ExtendPredictionForFlightPlan(vessel_guid);
  return m.Return(ToNewStatus(status));
}

void __cdecl principia__FlightPlanRenderedApsides(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const* const t_max,
    int const celestial_index,
    XYZ const sun_world_position,
    int const max_points,
    Iterator** const apoapsides,
    Iterator** const periapsides) {
  journal::Method<journal::FlightPlanRenderedApsides> m(
      {plugin,
       vessel_guid,
       t_max,
       celestial_index,
       sun_world_position,
       max_points},
      {apoapsides, periapsides});
  CHECK_NOTNULL(plugin);
  auto const& flight_plan =
      GetFlightPlan(*plugin, vessel_guid).GetAllSegments();
  DiscreteTrajectory<World> rendered_apoapsides;
  DiscreteTrajectory<World> rendered_periapsides;
  for (auto const& segment : flight_plan.segments()) {
    DiscreteTrajectory<World> segment_rendered_apoapsides;
    DiscreteTrajectory<World> segment_rendered_periapsides;
    plugin->ComputeAndRenderApsides(
        celestial_index,
        flight_plan,
        segment.begin(), segment.end(),
        t_max == nullptr ? InfiniteFuture : FromGameTime(*plugin, *t_max),
        FromXYZ<Position<World>>(sun_world_position),
        max_points,
        segment_rendered_apoapsides,
        segment_rendered_periapsides);
    rendered_apoapsides.Merge(std::move(segment_rendered_apoapsides));
    rendered_periapsides.Merge(std::move(segment_rendered_periapsides));
  }
  *apoapsides = new TypedIterator<DiscreteTrajectory<World>>(
      std::move(rendered_apoapsides),
      plugin);
  *periapsides = new TypedIterator<DiscreteTrajectory<World>>(
      std::move(rendered_periapsides),
      plugin);
  return m.Return();
}

void __cdecl principia__FlightPlanRenderedClosestApproaches(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position,
    int const max_points,
    Iterator** const closest_approaches) {
  journal::Method<journal::FlightPlanRenderedClosestApproaches> m(
      {plugin, vessel_guid, sun_world_position, max_points},
      {closest_approaches});
  CHECK_NOTNULL(plugin);
  auto const& flight_plan =
      GetFlightPlan(*plugin, vessel_guid).GetAllSegments();
  DiscreteTrajectory<World> rendered_closest_approaches;
  for (auto const& segment : flight_plan.segments()) {
    DiscreteTrajectory<World> segment_rendered_closest_approaches;
    plugin->ComputeAndRenderClosestApproaches(
        flight_plan,
        segment.begin(), segment.end(),
        FromXYZ<Position<World>>(sun_world_position),
        max_points,
        segment_rendered_closest_approaches);
    rendered_closest_approaches.Merge(
        std::move(segment_rendered_closest_approaches));
  }
  *closest_approaches = new TypedIterator<DiscreteTrajectory<World>>(
      std::move(rendered_closest_approaches),
      plugin);
  return m.Return();
}

void __cdecl principia__FlightPlanRenderedNodes(Plugin const* const plugin,
                                                char const* const vessel_guid,
                                                double const* const t_max,
                                                XYZ const sun_world_position,
                                                int const max_points,
                                                Iterator** const ascending,
                                                Iterator** const descending) {
  journal::Method<journal::FlightPlanRenderedNodes> m(
      {plugin, vessel_guid, t_max, sun_world_position, max_points},
      {ascending, descending});
  CHECK_NOTNULL(plugin);
  auto const& flight_plan =
      GetFlightPlan(*plugin, vessel_guid).GetAllSegments();
  std::vector<Renderer::Node> rendered_ascending;
  std::vector<Renderer::Node> rendered_descending;
  for (auto const& segment : flight_plan.segments()) {
    std::vector<Renderer::Node> segment_rendered_ascending;
    std::vector<Renderer::Node> segment_rendered_descending;
    plugin->ComputeAndRenderNodes(
        segment.begin(), segment.end(),
        t_max == nullptr ? InfiniteFuture : FromGameTime(*plugin, *t_max),
        FromXYZ<Position<World>>(sun_world_position),
        max_points,
        segment_rendered_ascending,
        segment_rendered_descending);
    std::move(segment_rendered_ascending.begin(),
              segment_rendered_ascending.end(),
              std::back_inserter(rendered_ascending));
    std::move(segment_rendered_descending.begin(),
              segment_rendered_descending.end(),
              std::back_inserter(rendered_descending));
  }
  *ascending = new TypedIterator<std::vector<Renderer::Node>>(
      std::move(rendered_ascending),
      plugin);
  *descending = new TypedIterator<std::vector<Renderer::Node>>(
      std::move(rendered_descending),
      plugin);
  return m.Return();
}

Iterator* __cdecl principia__FlightPlanRenderedSegment(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position,
    int const index) {
  journal::Method<journal::FlightPlanRenderedSegment> m({plugin,
                                                         vessel_guid,
                                                         sun_world_position,
                                                         index});
  CHECK_NOTNULL(plugin);

  // This might force a (partial) recomputation of the flight plan to avoid a
  // deadline, and a change of the anomalous status that will be noticed by the
  // flight planner.
  auto const segment =
      GetFlightPlan(*plugin, vessel_guid).GetSegmentAvoidingDeadlines(index);

  auto rendered_trajectory =
      plugin->renderer().RenderBarycentricTrajectoryInWorld(
          plugin->CurrentTime(),
          segment->begin(),
          segment->end(),
          FromXYZ<Position<World>>(sun_world_position),
          plugin->PlanetariumRotation());
  if (index % 2 == 1 && !rendered_trajectory.empty() &&
      rendered_trajectory.front().time != segment->front().time) {
    // TODO(egg): this is ugly; we should centralize rendering.
    // If this is a burn and we cannot render the beginning of the burn, we
    // render none of it, otherwise we try to render the Frenet trihedron at the
    // start and we fail.
    rendered_trajectory.clear();
  }
  return m.Return(new TypedIterator<DiscreteTrajectory<World>>(
      std::move(rendered_trajectory),
      plugin));
}

Status* __cdecl principia__FlightPlanReplace(Plugin const* const plugin,
                                             char const* const vessel_guid,
                                             Burn const& burn,
                                             int const index) {
  journal::Method<journal::FlightPlanReplace> m({plugin,
                                                 vessel_guid,
                                                 burn,
                                                 index});
  CHECK_NOTNULL(plugin);
  auto const status = GetFlightPlan(*plugin, vessel_guid)
                          .Replace(FromInterfaceBurn(*plugin, burn), index);
  plugin->ExtendPredictionForFlightPlan(vessel_guid);
  return m.Return(ToNewStatus(status));
}

void __cdecl principia__FlightPlanSelect(Plugin const* const plugin,
                                         char const* const vessel_guid,
                                         int const index) {
  journal::Method<journal::FlightPlanSelect> m({plugin, vessel_guid, index});
  CHECK_NOTNULL(plugin);
  plugin->GetVessel(vessel_guid)->SelectFlightPlan(index);
  return m.Return();
}

int __cdecl principia__FlightPlanSelected(Plugin const* const plugin,
                                          char const* const vessel_guid) {
  journal::Method<journal::FlightPlanSelected> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(plugin->GetVessel(vessel_guid)->selected_flight_plan_index());
}

Status* __cdecl principia__FlightPlanSetAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid,
    FlightPlanAdaptiveStepParameters const
        flight_plan_adaptive_step_parameters) {
  journal::Method<journal::FlightPlanSetAdaptiveStepParameters> m(
      {plugin, vessel_guid, flight_plan_adaptive_step_parameters});
  CHECK_NOTNULL(plugin);
  auto const parameters = FromFlightPlanAdaptiveStepParameters(
      flight_plan_adaptive_step_parameters);
  auto const status =
      GetFlightPlan(*plugin, vessel_guid)
          .SetAdaptiveStepParameters(parameters.first, parameters.second);
  plugin->ExtendPredictionForFlightPlan(vessel_guid);
  return m.Return(ToNewStatus(status));
}

Status* __cdecl principia__FlightPlanSetDesiredFinalTime(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const final_time) {
  journal::Method<journal::FlightPlanSetDesiredFinalTime> m({plugin,
                                                             vessel_guid,
                                                             final_time});
  CHECK_NOTNULL(plugin);
  auto const status =
      GetFlightPlan(*plugin, vessel_guid)
          .SetDesiredFinalTime(FromGameTime(*plugin, final_time));
  plugin->ExtendPredictionForFlightPlan(vessel_guid);
  return m.Return(ToNewStatus(status));
}

bool __cdecl principia__FlightPlanUpdateFromOptimization(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanUpdateFromOptimization> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(
      plugin->GetVessel(vessel_guid)->UpdateFlightPlanFromOptimization());
}

}  // namespace interface
}  // namespace principia
