
#include "ksp_plugin/interface.hpp"

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/vessel.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::not_null;
using geometry::Instant;
using ksp_plugin::Barycentric;
using ksp_plugin::FlightPlan;
using ksp_plugin::Navigation;
using ksp_plugin::NavigationManœuvre;
using ksp_plugin::Vessel;
using quantities::constants::StandardGravity;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using quantities::si::Tonne;

namespace interface {

namespace {

FlightPlan& GetFlightPlan(Plugin const* const plugin,
                          char const* const vessel_guid) {
  CHECK(CHECK_NOTNULL(plugin)->HasVessel(vessel_guid)) << vessel_guid;
  Vessel const& vessel = plugin->GetVessel(vessel_guid);
  CHECK(vessel.has_flight_plan()) << vessel_guid;
  return *vessel.flight_plan();
}

Burn GetBurn(NavigationManœuvre const& manœuvre) {
  return {manœuvre.thrust() / Kilo(Newton),
          manœuvre.specific_impulse() / (Second * StandardGravity),
          principia__GetNavigationFrameParameters(manœuvre.frame()),
          (manœuvre.initial_time() - Instant()) / Second,
          ToXYZ((manœuvre.Δv() / (Metre / Second)) *
                    manœuvre.direction().coordinates())};
}

ksp_plugin::Burn ToBurn(Plugin const* const plugin, Burn const& burn) {
  return {burn.thrust_in_kilonewtons * Kilo(Newton),
          burn.specific_impulse_in_seconds_g0 * Second * StandardGravity,
          base::check_not_null(std::unique_ptr<NavigationFrame>(
              principia__NewNavigationFrame(plugin, burn.frame))),
          Instant() + burn.initial_time * Second,
          Velocity<Frenet<Navigation>>(
              ToR3Element(burn.delta_v) * (Metre / Second))};
}

NavigationManoeuvre ToNavigationManoeuvre(NavigationManœuvre const& manœuvre) {
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
  result.direction = ToXYZ(manœuvre.direction().coordinates());
}

}  // namespace

bool principia__FlightPlanAppend(
    Plugin const* const plugin,
    char const* const vessel_guid,
    Burn const burn) {
  journal::Method<journal::FlightPlanAppend> m({plugin, vessel_guid, burn});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      Append(ToBurn(plugin, burn)));
}

NavigationManoeuvre principia__FlightPlanGetManoeuvre(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::FlightPlanGetManoeuvre> m({plugin,
                                                      vessel_guid,
                                                      index});
  return m.Return(ToNavigationManoeuvre(
             GetFlightPlan(plugin, vessel_guid).GetManœuvre(index)));
}

int principia__FlightPlanNumberOfManoeuvres(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfManoeuvres> m({plugin,
                                                            vessel_guid});
  return m.Return(GetFlightPlan(plugin, vessel_guid).number_of_manœuvres());
}

int principia__FlightPlanNumberOfSegments(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfSegments> m({plugin, vessel_guid});
  return m.Return(GetFlightPlan(plugin, vessel_guid).number_of_segments());
}

void principia__FlightPlanRemoveLast(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::FlightPlanRemoveLast> m({plugin, vessel_guid});
  GetFlightPlan(plugin, vessel_guid).RemoveLast();
  return m.Return();
}

LineAndIterator* principia__FlightPlanRenderedSegment(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::FlightPlanRenderedSegment> m({plugin,
                                                         vessel_guid,
                                                         index});
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(plugin, vessel_guid).GetSegment(index, &begin, &end);
  //Andthen?
}

bool principia__FlightPlanReplaceLast(
    Plugin const* const plugin,
    char const* const vessel_guid,
    Burn const burn) {
  journal::Method<journal::FlightPlanReplaceLast> m({plugin,
                                                     vessel_guid,
                                                     burn});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      ReplaceLast(ToBurn(plugin, burn)));
}

bool principia__FlightPlanSetFinalTime(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const final_time) {
  journal::Method<journal::FlightPlanSetFinalTime> m({plugin,
                                                      vessel_guid,
                                                      final_time});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      SetFinalTime(Instant() + final_time * Second));
}

void principia__FlightPlanSetTolerances(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const length_integration_tolerance,
    double const speed_integration_tolerance) {
  journal::Method<journal::FlightPlanSetTolerances>
      m({plugin,
         vessel_guid,
         length_integration_tolerance,
         speed_integration_tolerance});
  GetFlightPlan(plugin, vessel_guid).
      SetTolerances(length_integration_tolerance * Metre,
                    speed_integration_tolerance * (Metre / Second));
  return m.Return();
}

}  // namespace interface
}  // namespace principia
