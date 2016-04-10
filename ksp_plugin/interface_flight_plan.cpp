
#include "ksp_plugin/interface.hpp"

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::not_null;
using geometry::Instant;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using ksp_plugin::Barycentric;
using ksp_plugin::FlightPlan;
using ksp_plugin::Navigation;
using ksp_plugin::NavigationManœuvre;
using ksp_plugin::Positions;
using ksp_plugin::Vessel;
using ksp_plugin::World;
using ksp_plugin::WorldSun;
using physics::Ephemeris;
using quantities::constants::StandardGravity;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using quantities::si::Tonne;

namespace interface {

namespace {

Ephemeris<Barycentric>::AdaptiveStepParameters
FromInterfaceAdaptiveStepParameters(
    AdaptiveStepParameters const& adaptive_step_parameters) {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
      adaptive_step_parameters.max_steps,
      adaptive_step_parameters.length_integration_tolerance * Metre,
      adaptive_step_parameters.speed_integration_tolerance * (Metre / Second));
}

ksp_plugin::Burn FromInterfaceBurn(Plugin const* const plugin,
                                   Burn const& burn) {
  return {burn.thrust_in_kilonewtons * Kilo(Newton),
          burn.specific_impulse_in_seconds_g0 * Second * StandardGravity,
          NewNavigationFrame(plugin, burn.frame),
          Instant() + burn.initial_time * Second,
          Velocity<Frenet<Navigation>>(
              ToR3Element(burn.delta_v) * (Metre / Second))};
}

not_null<Vessel*> GetVessel(Plugin const* const plugin,
                            char const* const vessel_guid) {
  CHECK(CHECK_NOTNULL(plugin)->HasVessel(vessel_guid)) << vessel_guid;
  return plugin->GetVessel(vessel_guid);
}

FlightPlan& GetFlightPlan(Plugin const* const plugin,
                          char const* const vessel_guid) {
  Vessel const& vessel = *GetVessel(plugin, vessel_guid);
  CHECK(vessel.has_flight_plan()) << vessel_guid;
  return vessel.flight_plan();
}

Burn GetBurn(NavigationManœuvre const& manœuvre) {
  Velocity<Frenet<NavigationFrame>> const Δv =
      manœuvre.Δv() == Speed() ? Velocity<Frenet<NavigationFrame>>()
                               : manœuvre.Δv() * manœuvre.direction();

  // When building the parameters, make sure that the "optional" fields get a
  // deterministic default.
  NavigationFrameParameters parameters;
  parameters.centre_index = -1;
  parameters.primary_index = -1;
  parameters.secondary_index = -1;

  serialization::DynamicFrame message;
  manœuvre.frame()->WriteToMessage(&message);
  if (message.HasExtension(
          serialization::BarycentricRotatingDynamicFrame::
              barycentric_rotating_dynamic_frame)) {
    auto const& extension = message.GetExtension(
        serialization::BarycentricRotatingDynamicFrame::
            barycentric_rotating_dynamic_frame);
    parameters.extension = serialization::BarycentricRotatingDynamicFrame::
                               kBarycentricRotatingDynamicFrameFieldNumber;
    parameters.primary_index = extension.primary();
    parameters.secondary_index = extension.secondary();
  }
  if (message.HasExtension(
          serialization::BodyCentredNonRotatingDynamicFrame::
              body_centred_non_rotating_dynamic_frame)) {
    auto const& extension = message.GetExtension(
        serialization::BodyCentredNonRotatingDynamicFrame::
            body_centred_non_rotating_dynamic_frame);
    parameters.extension = serialization::BodyCentredNonRotatingDynamicFrame::
                               kBodyCentredNonRotatingDynamicFrameFieldNumber;
    parameters.centre_index = extension.centre();
  }

  return {manœuvre.thrust() / Kilo(Newton),
          manœuvre.specific_impulse() / (Second * StandardGravity),
          parameters,
          (manœuvre.initial_time() - Instant()) / Second,
          ToXYZ(Δv.coordinates() / (Metre / Second))};
}

AdaptiveStepParameters ToInterfaceAdaptiveStepParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  return {adaptive_step_parameters.max_steps(),
          adaptive_step_parameters.length_integration_tolerance() / Metre,
          adaptive_step_parameters.speed_integration_tolerance() /
              (Metre / Second)};
}

NavigationManoeuvre ToInterfaceNavigationManoeuvre(
    Plugin const* const plugin,
    NavigationManœuvre const& manœuvre) {
  OrthogonalMap<Barycentric, World> const barycentric_to_world =
      OrthogonalMap<WorldSun, World>::Identity() *
      plugin->BarycentricToWorldSun();
  NavigationManoeuvre result;
  result.burn = GetBurn(manœuvre);
  result.initial_mass_in_tonnes = manœuvre.initial_mass() / Tonne;
  result.final_mass_in_tonnes = manœuvre.final_mass() / Tonne;
  result.mass_flow = manœuvre.mass_flow() / (Kilogram / Second);
  result.duration = manœuvre.duration() / Second;
  result.final_time = (manœuvre.final_time() - Instant()) / Second;
  result.time_of_half_delta_v =
      (manœuvre.time_of_half_Δv() - Instant()) / Second;
  result.time_to_half_delta_v = manœuvre.time_to_half_Δv() / Second;
  Vector<double, Barycentric> const barycentric_inertial_direction =
      manœuvre.InertialDirection();
  Vector<double, World> const world_inertial_direction =
      barycentric_to_world(barycentric_inertial_direction);
  result.inertial_direction = ToXYZ(world_inertial_direction.coordinates());

  OrthogonalMap<Frenet<Navigation>, Barycentric> frenet_to_barycentric =
      manœuvre.FrenetFrame();
  Instant const current_time = plugin->CurrentTime();
  Instant const initial_time = manœuvre.initial_time();
  // TODO(egg): a separate |Frame| for plotted geometry (for now plotting goes
  // from |Barycentric| to itself, so it's easily forgotten), an utility for
  // plotting |Vector|s like the one for |Position|s.
  auto const plotting_frame = plugin->GetPlottingFrame();
  OrthogonalMap<Frenet<Navigation>, World> frenet_to_plotted_world =
      barycentric_to_world *
      plotting_frame->FromThisFrameAtTime(current_time).orthogonal_map() *
      plotting_frame->ToThisFrameAtTime(initial_time).orthogonal_map() *
      frenet_to_barycentric;
  result.tangent = ToXYZ(
      frenet_to_plotted_world(Vector<double, Frenet<Navigation>>({1, 0, 0}))
          .coordinates());
  result.normal = ToXYZ(
      frenet_to_plotted_world(Vector<double, Frenet<Navigation>>({0, 1, 0}))
          .coordinates());
  result.binormal = ToXYZ(
      frenet_to_plotted_world(Vector<double, Frenet<Navigation>>({0, 0, 1}))
          .coordinates());
  return result;
}

}  // namespace

bool principia__FlightPlanAppend(Plugin const* const plugin,
                                 char const* const vessel_guid,
                                 Burn const burn) {
  journal::Method<journal::FlightPlanAppend> m({plugin, vessel_guid, burn});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      Append(FromInterfaceBurn(plugin, burn)));
}

void principia__FlightPlanCreate(Plugin const* const plugin,
                                 char const* const vessel_guid,
                                 double const final_time,
                                 double const mass_in_tonnes) {
  journal::Method<journal::FlightPlanCreate> m({plugin,
                                                vessel_guid,
                                                final_time,
                                                mass_in_tonnes});
  CHECK_NOTNULL(plugin)->CreateFlightPlan(vessel_guid,
                                          Instant() + final_time * Second,
                                          mass_in_tonnes * Tonne);
  return m.Return();
}

void principia__FlightPlanDelete(Plugin const* const plugin,
                                 char const* const vessel_guid) {
  journal::Method<journal::FlightPlanDelete> m({plugin, vessel_guid});
  GetVessel(plugin, vessel_guid)->DeleteFlightPlan();
  return m.Return();
}

bool principia__FlightPlanExists(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanExists> m({plugin, vessel_guid});
  return m.Return(GetVessel(plugin, vessel_guid)->has_flight_plan());
}

AdaptiveStepParameters principia__FlightPlanGetAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetAdaptiveStepParameters> m(
      {plugin, vessel_guid});
  return m.Return(ToInterfaceAdaptiveStepParameters(
      GetFlightPlan(plugin, vessel_guid).adaptive_step_parameters()));
}

double principia__FlightPlanGetFinalTime(Plugin const* const plugin,
                                         char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetFinalTime> m({plugin, vessel_guid});
  return m.Return(
      (GetFlightPlan(plugin, vessel_guid).desired_final_time() - Instant()) /
      Second);
}

double principia__FlightPlanGetInitialTime(Plugin const* const plugin,
                                           char const* const vessel_guid) {
  journal::Method<journal::FlightPlanGetInitialTime> m({plugin, vessel_guid});
  return m.Return(
      (GetFlightPlan(plugin, vessel_guid).initial_time() - Instant()) / Second);
}

NavigationManoeuvre principia__FlightPlanGetManoeuvre(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::FlightPlanGetManoeuvre> m({plugin,
                                                      vessel_guid,
                                                      index});
  return m.Return(ToInterfaceNavigationManoeuvre(
                      plugin,
                      GetFlightPlan(plugin, vessel_guid).GetManœuvre(index)));
}

int principia__FlightPlanNumberOfManoeuvres(Plugin const* const plugin,
                                            char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfManoeuvres> m({plugin,
                                                            vessel_guid});
  return m.Return(GetFlightPlan(plugin, vessel_guid).number_of_manœuvres());
}

int principia__FlightPlanNumberOfSegments(Plugin const* const plugin,
                                          char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfSegments> m({plugin, vessel_guid});
  return m.Return(GetFlightPlan(plugin, vessel_guid).number_of_segments());
}

void principia__FlightPlanRemoveLast(Plugin const* const plugin,
                                     char const* const vessel_guid) {
  journal::Method<journal::FlightPlanRemoveLast> m({plugin, vessel_guid});
  GetFlightPlan(plugin, vessel_guid).RemoveLast();
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
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  CHECK_NOTNULL(plugin);
  GetFlightPlan(plugin, vessel_guid).GetAllSegments(&begin, &end);
  Position<World> q_sun =
      World::origin +
      Displacement<World>(ToR3Element(sun_world_position) * Metre);
  Positions<World> rendered_apoapsides;
  Positions<World> rendered_periapsides;
  plugin->ComputeAndRenderApsides(celestial_index,
                                  begin, end,
                                  q_sun,
                                  rendered_apoapsides,
                                  rendered_periapsides);
  *apoapsides = new TypedIterator<Positions<World>>(
      std::move(rendered_apoapsides));
  *periapsides = new TypedIterator<Positions<World>>(
      std::move(rendered_periapsides));
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
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(plugin, vessel_guid).GetSegment(index, &begin, &end);
  Positions<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedTrajectoryFromIterators(
          begin, end,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  return m.Return(new TypedIterator<Positions<World>>(
      std::move(rendered_trajectory)));
}

XYZSegment principia__FlightPlanRenderedSegmentEndpoints(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position,
    int const index) {
  journal::Method<journal::FlightPlanRenderedSegmentEndpoints> m(
      {plugin,
       vessel_guid,
       sun_world_position,
       index});
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(plugin, vessel_guid).GetSegment(index, &begin, &end);
  DiscreteTrajectory<Barycentric>::Iterator last = end;
  --last;
  Position<World> sun_position =
      World::origin +
      Displacement<World>(ToR3Element(sun_world_position) * Metre);
  Position<World> first_position =
      CHECK_NOTNULL(plugin)->PlotBarycentricPosition(
          begin.time(),
          begin.degrees_of_freedom().position(),
          sun_position);
  Position<World> last_position =
      CHECK_NOTNULL(plugin)->PlotBarycentricPosition(
          last.time(),
          last.degrees_of_freedom().position(),
          sun_position);
  XYZSegment result;
  result.begin = ToXYZ((first_position - World::origin).coordinates() / Metre);
  result.end = ToXYZ((last_position - World::origin).coordinates() / Metre);
  return m.Return(result);
}

bool principia__FlightPlanReplaceLast(Plugin const* const plugin,
                                      char const* const vessel_guid,
                                      Burn const burn) {
  journal::Method<journal::FlightPlanReplaceLast> m({plugin,
                                                     vessel_guid,
                                                     burn});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      ReplaceLast(FromInterfaceBurn(plugin, burn)));
}

bool principia__FlightPlanSetAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid,
    AdaptiveStepParameters const adaptive_step_parameters) {
  journal::Method<journal::FlightPlanSetAdaptiveStepParameters> m(
      {plugin, vessel_guid, adaptive_step_parameters});
  return m.Return(
      GetFlightPlan(plugin, vessel_guid).
          SetAdaptiveStepParameters(
              FromInterfaceAdaptiveStepParameters(adaptive_step_parameters)));
}

bool principia__FlightPlanSetFinalTime(Plugin const* const plugin,
                                       char const* const vessel_guid,
                                       double const final_time) {
  journal::Method<journal::FlightPlanSetFinalTime> m({plugin,
                                                      vessel_guid,
                                                      final_time});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      SetDesiredFinalTime(Instant() + final_time * Second));
}

}  // namespace interface
}  // namespace principia
