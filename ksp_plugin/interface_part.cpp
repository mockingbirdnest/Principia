
#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <chrono>

#include "geometry/grassmann.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/iterators.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using geometry::Bivector;
using quantities::Force;
using quantities::Torque;
using ksp_plugin::PartId;
using quantities::si::Kilo;
using quantities::si::Newton;
using quantities::si::Metre;
using quantities::si::Radian;

void __cdecl principia__PartApplyIntrinsicForce(
    Plugin* const plugin,
    PartId const part_id,
    XYZ const force_in_kilonewtons) {
  journal::Method<journal::PartApplyIntrinsicForce> m(
      {plugin, part_id, force_in_kilonewtons});
  CHECK_NOTNULL(plugin)->ApplyPartIntrinsicForce(
      part_id,
      Vector<Force, World>(FromXYZ(force_in_kilonewtons) * Kilo(Newton)));
  return m.Return();
}

void __cdecl principia__PartApplyIntrinsicForceAtPosition(
    Plugin* const plugin,
    PartId const part_id,
    XYZ const force_in_kilonewtons,
    XYZ const point_of_force_application,
    XYZ const part_position) {
  journal::Method<journal::PartApplyIntrinsicForceAtPosition> m(
      {plugin,
       part_id,
       force_in_kilonewtons,
       point_of_force_application,
       part_position});
  CHECK_NOTNULL(plugin)->ApplyPartIntrinsicForceAtPosition(
      part_id,
      Vector<Force, World>(FromXYZ(force_in_kilonewtons) * Kilo(Newton)),
      World::origin +
          Displacement<World>(FromXYZ(point_of_force_application) * Metre),
      World::origin + Displacement<World>(FromXYZ(part_position) * Metre));
  return m.Return();
}

void __cdecl principia__PartApplyIntrinsicTorque(
    Plugin* const plugin,
    PartId const part_id,
    XYZ const torque_in_kilonewton_metre) {
  journal::Method<journal::PartApplyIntrinsicTorque> m(
      {plugin, part_id, torque_in_kilonewton_metre});
  CHECK_NOTNULL(plugin)->ApplyPartIntrinsicTorque(
      part_id,
      Bivector<Torque, World>(FromXYZ(torque_in_kilonewton_metre) *
                              Kilo(Newton) * Metre * Radian));
  return m.Return();
}

QPRW __cdecl principia__PartGetActualDegreesOfFreedom(
    Plugin const* const plugin,
    PartId const part_id,
    Origin const origin) {
  journal::Method<journal::PartGetActualDegreesOfFreedom> m(
      {plugin, part_id, origin});
  CHECK_NOTNULL(plugin);
  RigidMotion<RigidPart, World> const part_motion = plugin->GetPartActualMotion(
      part_id,
      plugin->BarycentricToWorld(
          origin.reference_part_is_unmoving,
          origin.reference_part_id,
          origin.reference_part_is_at_origin
              ? std::nullopt
              : std::make_optional(FromXYZ<Position<World>>(
                    origin.main_body_centre_in_world))));
  DegreesOfFreedom<World> const part_dof =
      part_motion({RigidPart::origin, RigidPart::unmoving});
  Rotation<RigidPart, World> const part_orientation =
      part_motion.orthogonal_map().AsRotation();
  AngularVelocity<World> const part_angular_velocity =
      part_motion.Inverse().angular_velocity_of_to_frame();
  return m.Return(
      {ToQP(part_dof),
       ToWXYZ(part_orientation.quaternion()),
       ToXYZ(part_angular_velocity.coordinates() / (Radian / Second))});
}

bool __cdecl principia__PartIsTruthful(
    Plugin const* const plugin,
    uint32_t const part_id) {
  journal::Method<journal::PartIsTruthful> m({plugin, part_id});
  CHECK_NOTNULL(plugin);
  return m.Return(plugin->PartIsTruthful(part_id));
}

void __cdecl principia__PartSetApparentRigidMotion(
    Plugin* const plugin,
    PartId const part_id,
    QP const degrees_of_freedom,
    WXYZ const rotation,
    XYZ const angular_velocity) {
  journal::Method<journal::PartSetApparentRigidMotion> m(
      {plugin,
       part_id,
       degrees_of_freedom,
       rotation,
       angular_velocity});
  CHECK_NOTNULL(plugin);
  plugin->SetPartApparentRigidMotion(
      part_id,
      MakePartRigidMotion(degrees_of_freedom, rotation, angular_velocity));
  return m.Return();
}

}  // namespace interface
}  // namespace principia
