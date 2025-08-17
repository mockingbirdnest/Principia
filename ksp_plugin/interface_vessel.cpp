#include "ksp_plugin/interface.hpp"

#include <numeric>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/sign.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_sign;
using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_rigid_motion;
using namespace principia::quantities::_si;

XYZ __cdecl principia__VesselBinormal(Plugin const* const plugin,
                                      char const* const vessel_guid) {
  journal::Method<journal::VesselBinormal> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselBinormal(vessel_guid)));
}

// Calls `plugin->VesselFromParent` with the arguments given.
// `plugin` must not be null.  No transfer of ownership.
QP __cdecl principia__VesselFromParent(Plugin const* const plugin,
                                       int const parent_index,
                                       char const* const vessel_guid) {
  journal::Method<journal::VesselFromParent> m(
      {plugin, parent_index, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToQP(plugin->VesselFromParent(parent_index, vessel_guid)));
}

OrbitAnalysis* __cdecl principia__VesselGetAnalysis(
    Plugin* const plugin,
    char const* const vessel_guid,
    int const* const revolutions_per_cycle,
    int const* const days_per_cycle,
    int const ground_track_revolution) {
  journal::Method<journal::VesselGetAnalysis> m({plugin,
                                                 vessel_guid,
                                                 revolutions_per_cycle,
                                                 days_per_cycle,
                                                 ground_track_revolution});
  CHECK_NOTNULL(plugin);
  Vessel& vessel = *plugin->GetVessel(vessel_guid);
  vessel.RefreshOrbitAnalysis();
  not_null<OrbitAnalysis*> const analysis =
      NewOrbitAnalysis(vessel.orbit_analysis(),
                       *plugin,
                       revolutions_per_cycle,
                       days_per_cycle,
                       ground_track_revolution);
  analysis->progress_of_next_analysis = vessel.progress_of_orbit_analysis();
  return m.Return(analysis);
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

void __cdecl principia__VesselRequestAnalysis(Plugin* const plugin,
                                              char const* const vessel_guid,
                                              double const mission_duration) {
  journal::Method<journal::VesselRequestAnalysis> m(
      {plugin, vessel_guid, mission_duration});
  CHECK_NOTNULL(plugin);
  Vessel& vessel = *plugin->GetVessel(vessel_guid);
  plugin->ClearOrbitAnalysersOfVesselsOtherThan(vessel);
  vessel.RequestOrbitAnalysis(mission_duration * Second);
  return m.Return();
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

}  // namespace interface
}  // namespace principia
