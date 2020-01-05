#pragma once

#include "physics/closed_system.hpp"

namespace principia {
namespace physics {
namespace internal_closed_system {

using geometry::Displacement;
using geometry::OrthogonalMap;
using geometry::SymmetricProduct;
using geometry::Vector;
using geometry::Velocity;
using geometry::Wedge;
using physics::Anticommutator;
using quantities::Momentum;
using quantities::si::Radian;

template<typename InertialFrame, typename SystemFrame>
template<typename BodyFrame>
void ClosedSystem<InertialFrame, SystemFrame>::AddRigidBody(
    RigidMotion<BodyFrame, InertialFrame> motion,
    Mass mass,
    SymmetricBilinearForm<MomentOfInertia, BodyFrame> inertia_tensor) {
  DegreesOfFreedom<InertialFrame> degrees_of_freedom =
      motion({BodyFrame::origin, BodyFrame::unmoving});
  SymmetricBilinearForm<MomentOfInertia, InertialFrame>
      inertia_tensor_in_inertial_axes = motion.orthogonal_map()(inertia_tensor);
  centre_of_mass_.Add(degrees_of_freedom, mass);
  body_linear_motions_.emplace_back(degrees_of_freedom, mass);
  sum_of_inertia_tensors_ += inertia_tensor_in_inertial_axes;
  sum_of_intrinsic_angular_momenta_ +=
      Anticommutator(inertia_tensor_in_inertial_axes,
                     motion.Inverse().angular_velocity_of_to_frame());
}

template<typename InertialFrame, typename SystemFrame>
RigidMotion<SystemFrame, InertialFrame>
ClosedSystem<InertialFrame, SystemFrame>::linear_motion() const {
  DegreesOfFreedom<InertialFrame> centre_of_mass = centre_of_mass_.Get();
  return RigidMotion<SystemFrame, InertialFrame>(
      RigidTransformation<SystemFrame, InertialFrame>(
          SystemFrame::origin,
          centre_of_mass.position(),
          OrthogonalMap<SystemFrame, InertialFrame>::Identity()),
      InertialFrame::nonrotating,
      centre_of_mass.velocity());
}

template<typename InertialFrame, typename SystemFrame>
Mass ClosedSystem<InertialFrame, SystemFrame>::mass() const {
  return centre_of_mass_.weight();
}

template<typename InertialFrame, typename SystemFrame>
Bivector<AngularMomentum, SystemFrame>
ClosedSystem<InertialFrame, SystemFrame>::angular_momentum() const {
  RigidMotion<InertialFrame, SystemFrame> to_system_frame =
      linear_motion().Inverse();
  Bivector<AngularMomentum, SystemFrame> result =
      to_system_frame.orthogonal_map()(sum_of_intrinsic_angular_momenta_);
  for (auto const& [degrees_of_freedom, m] : body_linear_motions_) {
    DegreesOfFreedom<SystemFrame> const degrees_of_freedom_in_system_frame =
        to_system_frame(degrees_of_freedom);
    Displacement<SystemFrame> const r =
        degrees_of_freedom_in_system_frame.position() - SystemFrame::origin;
    Velocity<SystemFrame> const v =
        degrees_of_freedom_in_system_frame.velocity();
    Vector<Momentum, SystemFrame> const p = m * v;
    result += Wedge(r, p) * Radian;
  }
  return result;
}

template<typename InertialFrame, typename SystemFrame>
SymmetricBilinearForm<MomentOfInertia, SystemFrame>
ClosedSystem<InertialFrame, SystemFrame>::inertia_tensor() const {
  RigidMotion<InertialFrame, SystemFrame> to_system_frame =
      linear_motion().Inverse();
  SymmetricBilinearForm<MomentOfInertia, SystemFrame> result =
      to_system_frame.orthogonal_map()(sum_of_inertia_tensors_);
  for (auto const& [degrees_of_freedom, m] : body_linear_motions_) {
    DegreesOfFreedom<SystemFrame> const degrees_of_freedom_in_system_frame =
        to_system_frame(degrees_of_freedom);
    Displacement<SystemFrame> const r =
        degrees_of_freedom_in_system_frame.position() - SystemFrame::origin;
    result += m * SymmetricProduct(r, r);
  }
  return result;
}

}  // namespace internal_closed_system
}  // namespace physics
}  // namespace principia
