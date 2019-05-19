
#include "ksp_plugin/interface.hpp"

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/barycentric_rotating_dynamic_frame.hpp"
#include "physics/body_centred_body_direction_dynamic_frame.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using base::check_not_null;
using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Vector;
using geometry::Velocity;
using ksp_plugin::Barycentric;
using ksp_plugin::FlightPlan;
using ksp_plugin::Navigation;
using ksp_plugin::NavigationManœuvre;
using ksp_plugin::TypedIterator;
using ksp_plugin::Vessel;
using ksp_plugin::World;
using ksp_plugin::WorldSun;
using physics::BarycentricRotatingDynamicFrame;
using physics::BodyCentredBodyDirectionDynamicFrame;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::BodySurfaceDynamicFrame;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::Frenet;
using quantities::Speed;
using quantities::constants::StandardGravity;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using quantities::si::Tonne;

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

FlightPlan& GetFlightPlan(Plugin const& plugin,
                          char const* const vessel_guid) {
  Vessel const& vessel = *plugin.GetVessel(vessel_guid);
  CHECK(vessel.has_flight_plan()) << vessel_guid;
  return vessel.flight_plan();
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
    auto const* barycentric_rotating_dynamic_frame = dynamic_cast<
        BarycentricRotatingDynamicFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (barycentric_rotating_dynamic_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension =
          serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber;
      parameters.primary_index = plugin.CelestialIndexOfBody(
          *barycentric_rotating_dynamic_frame->primary());
      parameters.secondary_index = plugin.CelestialIndexOfBody(
          *barycentric_rotating_dynamic_frame->secondary());
    }
  }

  {
    auto const* body_centred_body_direction_dynamic_frame = dynamic_cast<
        BodyCentredBodyDirectionDynamicFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (body_centred_body_direction_dynamic_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension = serialization::
          BodyCentredBodyDirectionDynamicFrame::kExtensionFieldNumber;
      parameters.primary_index = plugin.CelestialIndexOfBody(
          *body_centred_body_direction_dynamic_frame->primary());
      parameters.secondary_index = plugin.CelestialIndexOfBody(
          *body_centred_body_direction_dynamic_frame->secondary());
    }
  }

  {
    auto const* body_centred_non_rotating_dynamic_frame = dynamic_cast<
        BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (body_centred_non_rotating_dynamic_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension = serialization::BodyCentredNonRotatingDynamicFrame::
          kExtensionFieldNumber;
      parameters.centre_index = plugin.CelestialIndexOfBody(
          *body_centred_non_rotating_dynamic_frame->centre());
    }
  }

  {
    auto const* body_surface_dynamic_frame = dynamic_cast<
        BodySurfaceDynamicFrame<Barycentric, Navigation> const*>(
            &*manœuvre.frame());
    if (body_surface_dynamic_frame != nullptr) {
      ++number_of_subclasses;
      parameters.extension =
          serialization::BodySurfaceDynamicFrame::kExtensionFieldNumber;
      parameters.centre_index =
          plugin.CelestialIndexOfBody(*body_surface_dynamic_frame->centre());
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

NavigationManoeuvre ToInterfaceNavigationManoeuvre(
    Plugin const& plugin,
    NavigationManœuvre const& manœuvre) {
  NavigationManoeuvre result;
  result.burn = GetBurn(plugin, manœuvre);
  result.initial_mass_in_tonnes = manœuvre.initial_mass() / Tonne;
  result.final_mass_in_tonnes = manœuvre.final_mass() / Tonne;
  result.mass_flow = manœuvre.mass_flow() / (Kilogram / Second);
  result.duration = manœuvre.duration() / Second;
  result.final_time = ToGameTime(plugin, manœuvre.final_time());
  result.time_of_half_delta_v = ToGameTime(plugin, manœuvre.time_of_half_Δv());
  result.time_to_half_delta_v = manœuvre.time_to_half_Δv() / Second;
  return result;
}

}  // namespace

bool principia__FlightPlanAppend(Plugin const* const plugin,
                                 char const* const vessel_guid,
                                 Burn const burn) {
  journal::Method<journal::FlightPlanAppend> m({plugin, vessel_guid, burn});
  CHECK_NOTNULL(plugin);

  // NOTE(phl): Preserving the previous semantics of FlightPlan.
  auto& flight_plan = GetFlightPlan(*plugin, vessel_guid);
  base::Status const status =
      flight_plan.Append(FromInterfaceBurn(*plugin, burn));
  if (status.error() == FlightPlan::singular ||
      status.error() == FlightPlan::does_not_fit) {
    return m.Return(false);
  } else if (status.ok() || flight_plan.number_of_anomalous_manœuvres() <= 1) {
    return m.Return(true);
  } else {
    flight_plan.RemoveLast();
    return m.Return(false);
  }
}

void principia__FlightPlanCreate(Plugin const* const plugin,
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

void principia__FlightPlanDelete(Plugin const* const plugin,
                                 char const* const vessel_guid) {
  journal::Method<journal::FlightPlanDelete> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  plugin->GetVessel(vessel_guid)->DeleteFlightPlan();
  return m.Return();
}

bool principia__FlightPlanExists(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanExists> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(plugin->GetVessel(vessel_guid)->has_flight_plan());
}

FlightPlanAdaptiveStepParameters
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

double principia__FlightPlanGetActualFinalTime(Plugin const* const plugin,
                                               char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetActualFinalTime> m(
      {plugin, vessel_guid});
  return m.Return(
      ToGameTime(*plugin,
                 GetFlightPlan(*plugin, vessel_guid).actual_final_time()));
}

double principia__FlightPlanGetDesiredFinalTime(Plugin const* const plugin,
                                                char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetDesiredFinalTime> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(
      ToGameTime(*plugin,
                 GetFlightPlan(*plugin, vessel_guid).desired_final_time()));
}

XYZ principia__FlightPlanGetGuidance(Plugin const* const plugin,
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
    result = plugin->renderer().FrenetToWorld(
                 *plugin->GetVessel(vessel_guid),
                 *manœuvre.frame(),
                 plugin->PlanetariumRotation())(manœuvre.direction());
  }
  return m.Return(ToXYZ(result));
}

double principia__FlightPlanGetInitialTime(Plugin const* const plugin,
                                           char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetInitialTime> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(
      ToGameTime(*plugin,
                 GetFlightPlan(*plugin, vessel_guid).initial_time()));
}

NavigationManoeuvre principia__FlightPlanGetManoeuvre(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::FlightPlanGetManoeuvre> m({plugin,
                                                      vessel_guid,
                                                      index});
  CHECK_NOTNULL(plugin);
  return m.Return(ToInterfaceNavigationManoeuvre(
                      *plugin,
                      GetFlightPlan(*plugin, vessel_guid).GetManœuvre(index)));
}

NavigationManoeuvreFrenetTrihedron
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

int principia__FlightPlanNumberOfManoeuvres(Plugin const* const plugin,
                                            char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfManoeuvres> m({plugin,
                                                            vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(GetFlightPlan(*plugin, vessel_guid).number_of_manœuvres());
}

int principia__FlightPlanNumberOfSegments(Plugin const* const plugin,
                                          char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfSegments> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(GetFlightPlan(*plugin, vessel_guid).number_of_segments());
}

void principia__FlightPlanRemoveLast(Plugin const* const plugin,
                                     char const* const vessel_guid) {
  journal::Method<journal::FlightPlanRemoveLast> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  GetFlightPlan(*plugin, vessel_guid).RemoveLast();
  return m.Return();
}

void principia__FlightPlanRenderedApsides(Plugin const* const plugin,
                                          char const* const vessel_guid,
                                          int const celestial_index,
                                          XYZ const sun_world_position,
                                          Iterator** const apoapsides,
                                          Iterator** const periapsides) {
  journal::Method<journal::FlightPlanRenderedApsides> m(
      {plugin, vessel_guid, celestial_index, sun_world_position},
      {apoapsides, periapsides});
  CHECK_NOTNULL(plugin);
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(*plugin, vessel_guid).GetAllSegments(begin, end);
  std::unique_ptr<DiscreteTrajectory<World>> rendered_apoapsides;
  std::unique_ptr<DiscreteTrajectory<World>> rendered_periapsides;
  plugin->ComputeAndRenderApsides(celestial_index,
                                  begin, end,
                                  FromXYZ<Position<World>>(sun_world_position),
                                  rendered_apoapsides,
                                  rendered_periapsides);
  *apoapsides = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_apoapsides)),
      plugin);
  *periapsides = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_periapsides)),
      plugin);
  return m.Return();
}

void principia__FlightPlanRenderedClosestApproaches(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position,
    Iterator** const closest_approaches) {
  journal::Method<journal::FlightPlanRenderedClosestApproaches> m(
      {plugin, vessel_guid, sun_world_position},
      {closest_approaches});
  CHECK_NOTNULL(plugin);
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(*plugin, vessel_guid).GetAllSegments(begin, end);
  std::unique_ptr<DiscreteTrajectory<World>> rendered_closest_approaches;
  plugin->ComputeAndRenderClosestApproaches(
      begin,
      end,
      FromXYZ<Position<World>>(sun_world_position),
      rendered_closest_approaches);
  *closest_approaches = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_closest_approaches)),
      plugin);
  return m.Return();
}

void principia__FlightPlanRenderedNodes(Plugin const* const plugin,
                                        char const* const vessel_guid,
                                        XYZ const sun_world_position,
                                        Iterator** const ascending,
                                        Iterator** const descending) {
  journal::Method<journal::FlightPlanRenderedNodes> m(
      {plugin, vessel_guid, sun_world_position},
      {ascending, descending});
  CHECK_NOTNULL(plugin);
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(*plugin, vessel_guid).GetAllSegments(begin, end);
  std::unique_ptr<DiscreteTrajectory<World>> rendered_ascending;
  std::unique_ptr<DiscreteTrajectory<World>> rendered_descending;
  plugin->ComputeAndRenderNodes(begin, end,
                                FromXYZ<Position<World>>(sun_world_position),
                                rendered_ascending,
                                rendered_descending);
  *ascending = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_ascending)),
      plugin);
  *descending = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_descending)),
      plugin);
  return m.Return();
}

Iterator* principia__FlightPlanRenderedSegment(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position,
    int const index) {
  journal::Method<journal::FlightPlanRenderedSegment> m({plugin,
                                                         vessel_guid,
                                                         sun_world_position,
                                                         index});
  CHECK_NOTNULL(plugin);
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(*plugin, vessel_guid).GetSegment(index, begin, end);
  auto rendered_trajectory =
      plugin->renderer().RenderBarycentricTrajectoryInWorld(
          plugin->CurrentTime(),
          begin,
          end,
          FromXYZ<Position<World>>(sun_world_position),
          plugin->PlanetariumRotation());
  if (index % 2 == 1 && !rendered_trajectory->Empty() &&
      rendered_trajectory->Begin().time() != begin.time()) {
    // TODO(egg): this is ugly; we should centralize rendering.
    // If this is a burn and we cannot render the beginning of the burn, we
    // render none of it, otherwise we try to render the Frenet trihedron at the
    // start and we fail.
    rendered_trajectory->ForgetAfter(astronomy::InfinitePast);
  }
  return m.Return(new TypedIterator<DiscreteTrajectory<World>>(
      std::move(rendered_trajectory),
      plugin));
}

bool principia__FlightPlanReplaceLast(Plugin const* const plugin,
                                      char const* const vessel_guid,
                                      Burn const burn) {
  journal::Method<journal::FlightPlanReplaceLast> m({plugin,
                                                     vessel_guid,
                                                     burn});
  CHECK_NOTNULL(plugin);

  // NOTE(phl): Preserving the previous semantics of FlightPlan.
  auto& flight_plan = GetFlightPlan(*plugin, vessel_guid);
  auto const manœuvre =
      flight_plan.GetManœuvre(flight_plan.number_of_manœuvres() - 1);
  base::Status const status =
      flight_plan.ReplaceLast(FromInterfaceBurn(*plugin, burn));
  if (status.error() == FlightPlan::singular ||
      status.error() == FlightPlan::does_not_fit) {
    return m.Return(false);
  } else if (status.ok() || flight_plan.number_of_anomalous_manœuvres() <= 1) {
    return m.Return(true);
  } else {
    flight_plan.Append(manœuvre.burn());
    return m.Return(false);
  }
}

bool principia__FlightPlanSetAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid,
    FlightPlanAdaptiveStepParameters const
        flight_plan_adaptive_step_parameters) {
  journal::Method<journal::FlightPlanSetAdaptiveStepParameters> m(
      {plugin, vessel_guid, flight_plan_adaptive_step_parameters});
  CHECK_NOTNULL(plugin);
  auto const parameters = FromFlightPlanAdaptiveStepParameters(
      flight_plan_adaptive_step_parameters);

  // NOTE(phl): Preserving the previous semantics of FlightPlan.
  auto& flight_plan = GetFlightPlan(*plugin, vessel_guid);
  auto const adaptive_step_parameters =
      flight_plan.adaptive_step_parameters();
  auto const generalized_adaptive_step_parameters =
      flight_plan.generalized_adaptive_step_parameters();
  base::Status const status =
      flight_plan.SetAdaptiveStepParameters(parameters.first,
                                            parameters.second);
  if (status.ok()) {
    return m.Return(true);
  } else {
    flight_plan.SetAdaptiveStepParameters(adaptive_step_parameters,
                                          generalized_adaptive_step_parameters);
    return m.Return(false);
  }
}

bool principia__FlightPlanSetDesiredFinalTime(Plugin const* const plugin,
                                              char const* const vessel_guid,
                                              double const final_time) {
  journal::Method<journal::FlightPlanSetDesiredFinalTime> m({plugin,
                                                             vessel_guid,
                                                             final_time});
  CHECK_NOTNULL(plugin);

  // NOTE(phl): Preserving the previous semantics of FlightPlan.
  auto& flight_plan = GetFlightPlan(*plugin, vessel_guid);
  base::Status const status =
      flight_plan.SetDesiredFinalTime(FromGameTime(*plugin, final_time));
  return m.Return(status.error() != FlightPlan::bad_desired_final_time);
}

}  // namespace interface
}  // namespace principia
