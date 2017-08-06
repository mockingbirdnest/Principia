
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

void principia__VesselSetPredictionAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid,
    AdaptiveStepParameters const adaptive_step_parameters) {
  journal::Method<journal::VesselSetPredictionAdaptiveStepParameters> m(
      {plugin, vessel_guid, adaptive_step_parameters});
  CHECK_NOTNULL(plugin);
  plugin->GetVessel(vessel_guid)
      ->set_prediction_adaptive_step_parameters(
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
