
#include "ksp_plugin/interface.hpp"

#include <numeric>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::OrthogonalMap;
using geometry::RigidTransformation;
using geometry::Sign;
using geometry::Vector;
using geometry::Velocity;
using ksp_plugin::AliceSun;
using ksp_plugin::Apparent;
using physics::DegreesOfFreedom;
using physics::RelativeDegreesOfFreedom;
using physics::RigidMotion;
using quantities::Force;
using quantities::si::Kilo;
using quantities::si::Newton;
using quantities::si::Tonne;

XYZ __cdecl principia__VesselBinormal(Plugin const* const plugin,
                                      char const* const vessel_guid) {
  journal::Method<journal::VesselBinormal> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselBinormal(vessel_guid)));
}

// Calls |plugin->VesselFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
QP __cdecl principia__VesselFromParent(Plugin const* const plugin,
                                       int const parent_index,
                                       char const* const vessel_guid) {
  journal::Method<journal::VesselFromParent> m(
      {plugin, parent_index, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToQP(plugin->VesselFromParent(parent_index, vessel_guid)));
}

AdaptiveStepParameters __cdecl
principia__VesselGetPredictionAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::VesselGetPredictionAdaptiveStepParameters> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToAdaptiveStepParameters(
      plugin->GetVessel(vessel_guid)->prediction_adaptive_step_parameters()));
}

XYZ __cdecl principia__VesselNormal(Plugin const* const plugin,
                                    char const* const vessel_guid) {
  journal::Method<journal::VesselNormal> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselNormal(vessel_guid)));
}

OrbitAnalysis* __cdecl principia__VesselRefreshAnalysis(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const mission_duration,
    int const* const revolutions_per_cycle,
    int const* const days_per_cycle,
    int const ground_track_revolution) {
  journal::Method<journal::VesselRefreshAnalysis> m({plugin,
                                                     vessel_guid,
                                                     mission_duration,
                                                     revolutions_per_cycle,
                                                     days_per_cycle,
                                                     ground_track_revolution});
  CHECK_NOTNULL(plugin);
  CHECK_EQ(revolutions_per_cycle == nullptr, days_per_cycle == nullptr);
  bool const has_nominal_recurrence = revolutions_per_cycle != nullptr;
  if (has_nominal_recurrence) {
    CHECK_GT(*revolutions_per_cycle, 0);
    CHECK_NE(*days_per_cycle, 0);
  }
  Vessel& vessel = *plugin->GetVessel(vessel_guid);
  vessel.RefreshOrbitAnalysis(mission_duration * Second);
  OrbitAnalysis* const analysis = new OrbitAnalysis{};
  analysis->progress_of_next_analysis = vessel.progress_of_orbit_analysis();
  if (vessel.orbit_analysis() != nullptr) {
    analysis->primary_index = new int(
        plugin->CelestialIndexOfBody(vessel.orbit_analysis()->primary()));
    analysis->mission_duration =
        vessel.orbit_analysis()->mission_duration() / Second;
    if (vessel.orbit_analysis()->elements().has_value()) {
      auto const& elements = *vessel.orbit_analysis()->elements();
      analysis->elements = new OrbitalElements{
          .sidereal_period = elements.sidereal_period() / Second,
          .nodal_period = elements.nodal_period() / Second,
          .anomalistic_period = elements.anomalistic_period() / Second,
          .nodal_precession = elements.nodal_precession() / (Radian / Second),
          .mean_semimajor_axis =
          ToInterval(elements.mean_semimajor_axis_interval()),
          .mean_eccentricity =
          ToInterval(elements.mean_eccentricity_interval()),
          .mean_inclination = ToInterval(elements.mean_inclination_interval()),
          .mean_longitude_of_ascending_nodes =
          ToInterval(elements.mean_longitude_of_ascending_node_interval()),
          .mean_argument_of_periapsis =
          ToInterval(elements.mean_argument_of_periapsis_interval()),
      };
    }
    if (has_nominal_recurrence) {
      int const Cᴛₒ =
          Sign(vessel.orbit_analysis()->primary().angular_frequency()) *
          std::abs(*days_per_cycle);
      int const νₒ =
          std::nearbyint(static_cast<double>(*revolutions_per_cycle) / Cᴛₒ);
      int const Dᴛₒ = *revolutions_per_cycle - νₒ * Cᴛₒ;
      int const gcd = std::gcd(Dᴛₒ, Cᴛₒ);
      vessel.orbit_analysis()->SetRecurrence({νₒ, Dᴛₒ / gcd, Cᴛₒ / gcd});
    } else {
      vessel.orbit_analysis()->ResetRecurrence();
    }
    if (vessel.orbit_analysis()->recurrence().has_value()) {
      auto const& recurrence = *vessel.orbit_analysis()->recurrence();
      analysis->recurrence = new OrbitRecurrence{
          .nuo = recurrence.νₒ(),
          .dto = recurrence.Dᴛₒ(),
          .cto = recurrence.Cᴛₒ(),
          .number_of_revolutions =
          recurrence.number_of_revolutions(),
          .equatorial_shift =
          recurrence.equatorial_shift() / Radian,
          .base_interval = recurrence.base_interval() / Radian,
          .grid_interval = recurrence.grid_interval() / Radian,
          .subcycle = recurrence.subcycle(),
      };
    }
    if (vessel.orbit_analysis()->ground_track().has_value() &&
        vessel.orbit_analysis()->equatorial_crossings().has_value()) {
      auto const& equatorial_crossings =
          *vessel.orbit_analysis()->equatorial_crossings();
      analysis->ground_track = new OrbitGroundTrack{
          .equatorial_crossings = {
              .longitudes_reduced_to_ascending_pass =
              ToInterval(equatorial_crossings.longitudes_reduced_to_pass(
                  2 * ground_track_revolution - 1)),
              .longitudes_reduced_to_descending_pass =
              ToInterval(equatorial_crossings.longitudes_reduced_to_pass(
                  2 * ground_track_revolution)),
          },
      };
    }
  }

  return m.Return(analysis);
}

void __cdecl principia__VesselSetPredictionAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid,
    AdaptiveStepParameters const adaptive_step_parameters) {
  journal::Method<journal::VesselSetPredictionAdaptiveStepParameters> m(
      {plugin, vessel_guid, adaptive_step_parameters});
  CHECK_NOTNULL(plugin);
  plugin->SetPredictionAdaptiveStepParameters(
      vessel_guid, FromAdaptiveStepParameters(adaptive_step_parameters));
  return m.Return();
}

XYZ __cdecl principia__VesselTangent(Plugin const* const plugin,
                                     char const* const vessel_guid) {
  journal::Method<journal::VesselTangent> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselTangent(vessel_guid)));
}

XYZ __cdecl principia__VesselVelocity(Plugin const* const plugin,
                                      char const* const vessel_guid) {
  journal::Method<journal::VesselVelocity> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselVelocity(vessel_guid)));
}

void __cdecl principia__SetAngularMomentumConservation(
    bool const conserve_angular_momentum) {
  journal::Method<journal::SetAngularMomentumConservation> m(
      {conserve_angular_momentum});
  ksp_plugin::PileUp::conserve_angular_momentum = conserve_angular_momentum;
  return m.Return();
}

char const* __cdecl principia__VesselGetPileUpTrace(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::VesselGetPileUpTrace> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  char const* trace;
  plugin->GetVessel(vessel_guid)->ForSomePart([&trace](ksp_plugin::Part& part) {
    trace = part.containing_pile_up()->trace.data();
  });
  return m.Return(trace);
}

}  // namespace interface
}  // namespace principia
