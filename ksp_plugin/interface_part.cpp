#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>

#include "absl/log/check.h"
#include "absl/log/die_if_null.h"
#include "absl/log/log.h"
#include "base/algebra.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/part.hpp"
#include "physics/tensors.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::base::_algebra;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_identification;
using namespace principia::ksp_plugin::_part;
using namespace principia::physics::_tensors;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

namespace {

InertiaTensor<RigidPart> FromMomentsOfInertia(
    XYZ const& moments_of_inertia_in_tonnes,
    WXYZ const& principal_axes_rotation) {
  static constexpr MomentOfInertia zero;
  auto const moments_of_inertia = FromXYZ<R3Element<MomentOfInertia>>(
      {.x = moments_of_inertia_in_tonnes.x,
       .y = moments_of_inertia_in_tonnes.y,
       .z = moments_of_inertia_in_tonnes.z});
  InertiaTensor<PartPrincipalAxes> const inertia_tensor_in_principal_axes(
      R3x3Matrix<MomentOfInertia>({moments_of_inertia.x, zero, zero},
                                  {zero, moments_of_inertia.y, zero},
                                  {zero, zero, moments_of_inertia.z}));

  Rotation<PartPrincipalAxes, RigidPart> const principal_axes_to_part(
      FromWXYZ(principal_axes_rotation));

  return principal_axes_to_part(inertia_tensor_in_principal_axes);
}

}  // namespace

XYZ __cdecl principia__PartAngularMomentumFromAngularVelocity(
    XYZ world_angular_velocity,
    XYZ moments_of_inertia_in_tonnes,
    WXYZ principal_axes_rotation,
    WXYZ part_rotation) {
  journal::Method<journal::PartAngularMomentumFromAngularVelocity> m(
      {world_angular_velocity,
       moments_of_inertia_in_tonnes,
       principal_axes_rotation,
       part_rotation});
  auto const angular_velocity =
      FromXYZ<AngularVelocity<World>>(world_angular_velocity);
  Rotation<RigidPart, World> const part_to_world(FromWXYZ(part_rotation));
  InertiaTensor<World> const inertia_tensor =
      part_to_world(FromMomentsOfInertia(moments_of_inertia_in_tonnes,
                                         principal_axes_rotation));

  Bivector<AngularMomentum, World> const angular_momentum =
      inertia_tensor * angular_velocity;

  return m.Return(ToXYZ(angular_momentum));
}

void __cdecl principia__PartApplyIntrinsicForce(
    Plugin* const plugin,
    PartId const part_id,
    XYZ const force_in_kilonewtons) {
  journal::Method<journal::PartApplyIntrinsicForce> m(
      {plugin, part_id, force_in_kilonewtons});
  ABSL_DIE_IF_NULL(plugin)->ApplyPartIntrinsicForce(
      part_id,
      Vector<Force, World>(FromXYZ(force_in_kilonewtons) * Kilo(Newton)));
  return m.Return();
}

void __cdecl principia__PartApplyIntrinsicForceAtPosition(
    Plugin* const plugin,
    PartId const part_id,
    XYZ const force_in_kilonewtons,
    XYZ const lever_arm) {
  journal::Method<journal::PartApplyIntrinsicForceAtPosition> m(
      {plugin,
       part_id,
       force_in_kilonewtons,
       lever_arm});
  ABSL_DIE_IF_NULL(plugin)->ApplyPartIntrinsicForceAtPosition(
      part_id,
      Vector<Force, World>(FromXYZ(force_in_kilonewtons) * Kilo(Newton)),
      Displacement<World>(FromXYZ(lever_arm) * Metre));
  return m.Return();
}

void __cdecl principia__PartApplyIntrinsicTorque(
    Plugin* const plugin,
    PartId const part_id,
    XYZ const torque_in_kilonewton_metre) {
  journal::Method<journal::PartApplyIntrinsicTorque> m(
      {plugin, part_id, torque_in_kilonewton_metre});
  ABSL_DIE_IF_NULL(plugin)->ApplyPartIntrinsicTorque(
      part_id,
      Bivector<Torque, World>(FromXYZ(torque_in_kilonewton_metre) *
                              Kilo(Newton) * Metre * Radian));
  return m.Return();
}

XYZ __cdecl principia__PartDragTorqueFromAngularVelocity(
    double const angular_drag_per_second,
    double const delta_t,
    XYZ const world_angular_velocity,
    XYZ const moments_of_inertia_in_tonnes,
    WXYZ const principal_axes_rotation,
    WXYZ const part_rotation) {
  journal::Method<journal::PartDragTorqueFromAngularVelocity> m(
      {angular_drag_per_second,
       delta_t,
       world_angular_velocity,
       moments_of_inertia_in_tonnes,
       principal_axes_rotation,
       part_rotation});
  Time const Δt = delta_t * Second;
  // The type of the angular drag coefficient is given by its use in PhysX, see
  // https://github.com/NVIDIAGameWorks/PhysX-3.4/blob/5e42a5f112351a223c19c17bb331e6c55037b8eb/PhysX_3.4/Source/LowLevelDynamics/src/DyBodyCoreIntegrator.h#L72-L75.
  Inverse<Time> const angular_drag = angular_drag_per_second / Second;

  Rotation<RigidPart, World> const part_to_world(FromWXYZ(part_rotation));
  Rotation<World, RigidPart> const world_to_part = part_to_world.Inverse();

  AngularVelocity<RigidPart> const angular_velocity =
      world_to_part(FromXYZ<AngularVelocity<World>>(world_angular_velocity));
  InertiaTensor<RigidPart> const inertia_tensor = FromMomentsOfInertia(
      moments_of_inertia_in_tonnes, principal_axes_rotation);

  return m.Return(ToXYZ(part_to_world(Part::DragTorqueFromAngularVelocity(
      angular_drag, Δt, angular_velocity, inertia_tensor))));
}

QPRW __cdecl principia__PartGetActualRigidMotion(
    Plugin const* const plugin,
    PartId const part_id,
    Origin const origin) {
  journal::Method<journal::PartGetActualRigidMotion> m(
      {plugin, part_id, origin});
  CHECK(plugin != nullptr);
  RigidMotion<EccentricPart, World> const part_motion =
      plugin->GetPartActualMotion(
          part_id,
          plugin->BarycentricToWorld(
              origin.reference_part_is_unmoving,
              origin.reference_part_id,
              origin.reference_part_is_at_origin
                  ? std::nullopt
                  : std::make_optional(FromXYZ<Position<World>>(
                        origin.main_body_centre_in_world))));
  DegreesOfFreedom<World> const part_dof =
      part_motion({EccentricPart::origin, EccentricPart::unmoving});
  Rotation<EccentricPart, World> const part_orientation =
      part_motion.orthogonal_map().AsRotation();
  AngularVelocity<World> const part_angular_velocity =
      part_motion.angular_velocity_of<EccentricPart>();
  return m.Return(
      {ToQP(part_dof),
       ToWXYZ(part_orientation.quaternion()),
       ToXYZ(part_angular_velocity.coordinates() / (Radian / Second))});
}

bool __cdecl principia__PartIsTruthful(
    Plugin const* const plugin,
    uint32_t const part_id) {
  journal::Method<journal::PartIsTruthful> m({plugin, part_id});
  CHECK(plugin != nullptr);
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
  CHECK(plugin != nullptr);
  plugin->SetPartApparentRigidMotion(
      part_id,
      MakePartApparentRigidMotion(
          degrees_of_freedom, rotation, angular_velocity));
  return m.Return();
}

}  // namespace interface
}  // namespace principia
