
#include "ksp_plugin/interface.hpp"

#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::OrthogonalMap;
using geometry::RigidTransformation;
using geometry::Vector;
using geometry::Velocity;
using ksp_plugin::AliceSun;
using ksp_plugin::ApparentBubble;
using physics::DegreesOfFreedom;
using physics::RelativeDegreesOfFreedom;
using physics::RigidMotion;
using quantities::Force;
using quantities::si::Kilo;
using quantities::si::Newton;
using quantities::si::Tonne;

XYZ principia__VesselBinormal(Plugin const* const plugin,
                              char const* const vessel_guid) {
  journal::Method<journal::VesselBinormal> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselBinormal(vessel_guid)));
}

// Calls |plugin->VesselFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
QP principia__VesselFromParent(Plugin const* const plugin,
                               int const parent_index,
                               char const* const vessel_guid) {
  journal::Method<journal::VesselFromParent> m(
      {plugin, parent_index, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToQP(plugin->VesselFromParent(parent_index, vessel_guid)));
}

AdaptiveStepParameters principia__VesselGetPredictionAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::VesselGetPredictionAdaptiveStepParameters> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToAdaptiveStepParameters(
      plugin->GetVessel(vessel_guid)->prediction_adaptive_step_parameters()));
}

XYZ principia__VesselNormal(Plugin const* const plugin,
                            char const* const vessel_guid) {
  journal::Method<journal::VesselNormal> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselNormal(vessel_guid)));
}

OrbitAnalysis principia__VesselRefreshAnalysis(Plugin const* const plugin,
                                               char const* const vessel_guid,
                                               int const primary_index,
                                               double const mission_duration) {
  journal::Method<journal::VesselRefreshAnalysis> m({plugin, vessel_guid, primary_index});
  CHECK_NOTNULL(plugin);
  Vessel& vessel = *plugin->GetVessel(vessel_guid);
  vessel.RefreshOrbitAnalysis(plugin->GetCelestial(primary_index).body(),
                              mission_duration * Second);
  OrbitAnalysis analysis{};
  analysis.progress_percentage = vessel.orbit_analysis_percentage();
  if (vessel.orbit_analysis().has_value()) {
    analysis.primary_index =
        plugin->CelestialIndexOfBody(*vessel.orbit_analysis()->primary);
    analysis.mission_duration =
        vessel.orbit_analysis()->mission_duration / Second;
    if (vessel.orbit_analysis()->elements.has_value()) {
      auto const& elements = *vessel.orbit_analysis()->elements;
      analysis.elements.anomalistic_period =
          elements.anomalistic_period() / Second;
      analysis.elements.nodal_period = elements.nodal_period() / Second;
      analysis.elements.sidereal_period = elements.sidereal_period() / Second;
      analysis.elements.mean_argument_of_periapsis =
          ToInterval(elements.mean_argument_of_periapsis_interval());
      analysis.elements.mean_eccentricity =
          ToInterval(elements.mean_eccentricity_interval());
      analysis.elements.mean_inclination =
          ToInterval(elements.mean_inclination_interval());
      analysis.elements.mean_longitude_of_ascending_nodes =
          ToInterval(elements.mean_longitude_of_ascending_node_interval());
      analysis.elements.mean_semimajor_axis =
          ToInterval(elements.mean_semimajor_axis_interval());
    }
    if (vessel.orbit_analysis()->recurrence.has_value()) {
      auto const& recurrence = *vessel.orbit_analysis()->recurrence;
      analysis.recurrence.nuo = recurrence.νₒ();
      analysis.recurrence.dto = recurrence.Dᴛₒ();
      analysis.recurrence.cto = recurrence.Cᴛₒ();
      analysis.recurrence.number_of_revolutions =
          recurrence.number_of_revolutions();
      analysis.recurrence.subcycle = recurrence.subcycle();
      analysis.recurrence.equatorial_shift =
          recurrence.equatorial_shift() / Radian;
      analysis.recurrence.base_interval = recurrence.base_interval() / Radian;
      analysis.recurrence.grid_interval = recurrence.grid_interval() / Radian;
    }
    if (vessel.orbit_analysis()->ground_track.has_value()) {
      auto const& ground_track = *vessel.orbit_analysis()->ground_track;
      if (auto const& longitudes =
              ground_track
                  .reduced_longitudes_of_equator_crossings_of_ascending_passes();
          longitudes.has_value()) {
        analysis.ground_track
            .reduced_longitudes_of_equator_crossings_of_ascending_passes =
            ToInterval(*longitudes);
      }
      if (auto const& longitudes =
              ground_track
                  .reduced_longitudes_of_equator_crossings_of_descending_passes();
          longitudes.has_value()) {
        analysis.ground_track
            .reduced_longitudes_of_equator_crossings_of_descending_passes =
            ToInterval(*longitudes);
      }
    }
  }

  return m.Return(analysis);
}

void principia__VesselSetPredictionAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid,
    AdaptiveStepParameters const adaptive_step_parameters) {
  journal::Method<journal::VesselSetPredictionAdaptiveStepParameters> m(
      {plugin, vessel_guid, adaptive_step_parameters});
  CHECK_NOTNULL(plugin);
  plugin->SetPredictionAdaptiveStepParameters(
      vessel_guid,
      FromAdaptiveStepParameters(adaptive_step_parameters));
  return m.Return();
}

XYZ principia__VesselTangent(Plugin const* const plugin,
                             char const* const vessel_guid) {
  journal::Method<journal::VesselTangent> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselTangent(vessel_guid)));
}

XYZ principia__VesselVelocity(Plugin const* const plugin,
                              char const* const vessel_guid) {
  journal::Method<journal::VesselVelocity> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselVelocity(vessel_guid)));
}

}  // namespace interface
}  // namespace principia
