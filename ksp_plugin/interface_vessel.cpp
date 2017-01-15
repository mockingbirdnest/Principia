
#include "ksp_plugin/interface.hpp"

#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "geometry/grassmann.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using geometry::Vector;
using quantities::Force;
using quantities::si::Kilo;
using quantities::si::Newton;
using quantities::si::Tonne;

XYZ principia__VesselBinormal(Plugin const* const plugin,
                              char const* const vessel_guid) {
  journal::Method<journal::VesselBinormal> m({plugin, vessel_guid});
  return m.Return(
      ToXYZ(CHECK_NOTNULL(plugin)->VesselBinormal(vessel_guid).coordinates()));
}

void principia__VesselClearIntrinsicForce(Plugin const* const plugin,
                                          char const* const vessel_guid) {
  journal::Method<journal::VesselClearIntrinsicForce> m({plugin, vessel_guid});
  GetVessel(*plugin, vessel_guid)->clear_intrinsic_force();
  return m.Return();
}

void principia__VesselClearMass(Plugin const* const plugin,
                                char const* const vessel_guid) {
  journal::Method<journal::VesselClearMass> m({plugin, vessel_guid});
  GetVessel(*plugin, vessel_guid)->clear_mass();
  return m.Return();
}

AdaptiveStepParameters principia__VesselGetPredictionAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid) {
  journal::Method<journal::VesselGetPredictionAdaptiveStepParameters> m(
      {plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToAdaptiveStepParameters(
      GetVessel(*plugin, vessel_guid)->prediction_adaptive_step_parameters()));
}

void principia__VesselIncrementIntrinsicForce(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const intrinsic_force_in_kilonewtons) {
  journal::Method<journal::VesselIncrementIntrinsicForce> m(
      {plugin, vessel_guid, intrinsic_force_in_kilonewtons});
  // TODO(phl): this is profoundly incorrect! the vector is given in World, not
  // Barycentric (KSP doesn't know about Barycentric).
  GetVessel(*plugin, vessel_guid)
      ->increment_intrinsic_force(Vector<Force, Barycentric>(
          FromXYZ(intrinsic_force_in_kilonewtons) * Kilo(Newton)));
  return m.Return();
}

void principia__VesselIncrementMass(Plugin const* const plugin,
                                    char const* const vessel_guid,
                                    double const mass_in_tonnes) {
  journal::Method<journal::VesselIncrementMass> m(
      {plugin, vessel_guid, mass_in_tonnes});
  GetVessel(*plugin, vessel_guid)->increment_mass(mass_in_tonnes * Tonne);
  return m.Return();
}

XYZ principia__VesselNormal(Plugin const* const plugin,
                            char const* const vessel_guid) {
  journal::Method<journal::VesselNormal> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselNormal(vessel_guid).coordinates()));
}

void principia__VesselSetPredictionAdaptiveStepParameters(
    Plugin const* const plugin,
    char const* const vessel_guid,
    AdaptiveStepParameters const adaptive_step_parameters) {
  journal::Method<journal::VesselSetPredictionAdaptiveStepParameters> m(
      {plugin, vessel_guid, adaptive_step_parameters});
  CHECK_NOTNULL(plugin);
  GetVessel(*plugin, vessel_guid)
      ->set_prediction_adaptive_step_parameters(
          FromAdaptiveStepParameters(adaptive_step_parameters));
  return m.Return();
}

XYZ principia__VesselTangent(Plugin const* const plugin,
                             char const* const vessel_guid) {
  journal::Method<journal::VesselTangent> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselTangent(vessel_guid).coordinates()));
}

XYZ principia__VesselVelocity(Plugin const* const plugin,
                              char const* const vessel_guid) {
  journal::Method<journal::VesselVelocity> m({plugin, vessel_guid});
  CHECK_NOTNULL(plugin);
  return m.Return(ToXYZ(plugin->VesselVelocity(vessel_guid).coordinates() /
                        (Metre / Second)));
}

}  // namespace interface
}  // namespace principia
