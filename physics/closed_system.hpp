#pragma once

#include <optional>
#include <vector>

#include "physics/inertia_tensor.hpp"
#include "physics/rigid_motion.hpp"

namespace principia {
namespace physics {
namespace internal_closed_system {

using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Frame;
using geometry::OrthogonalMap;
using geometry::SymmetricBilinearForm;
using quantities::AngularMomentum;
using quantities::Mass;
using quantities::MomentOfInertia;

// Computes the instantaneous overall properties of a mechanical system.
// Effectively an extension of
// |BarycentreCalculator<DegreesOfFreedom<InertialFrame>, Mass>| that also
// handles rotational motion.
template<typename InertialFrame, typename SystemFrame>
class ClosedSystem {
 public:
  ClosedSystem() = default;

  template<typename BodyFrame>
  void AddRigidBody(
      RigidMotion<BodyFrame, InertialFrame> motion,
      Mass mass,
      SymmetricBilinearForm<MomentOfInertia, BodyFrame> inertia_tensor);

  // The motion of a non-rotating frame whose origin is the centre of mass of
  // the system, and whose axes are those of |InertialFrame|.
  RigidMotion<SystemFrame, InertialFrame> linear_motion() const;
  // The total mass.
  Mass mass() const;
  // The total angular momentum of the system, with respect to the origin of
  // |SystemFrame|, i.e., the origin of the system.
  Bivector<AngularMomentum, SystemFrame> angular_momentum() const;
  // The moments of inertia of the system, with respect to the origin of
  // |SystemFrame|, i.e., the origin of the system.
  SymmetricBilinearForm<MomentOfInertia, SystemFrame> inertia_tensor() const;

 private:
  BarycentreCalculator<DegreesOfFreedom<InertialFrame>, Mass> centre_of_mass_;
  std::vector<std::pair<DegreesOfFreedom<InertialFrame>, Mass>>
      body_linear_motions_;
  // The sum of the intrinsic angular momenta of the bodies, i.e., of the
  // angular momenta with respect to their individual centres of mass.  This is
  // not the total angular momentum, to which the linear motions of the bodies
  // with respect to the centre of mass of the system also contribute.
  Bivector<AngularMomentum, InertialFrame> sum_of_intrinsic_angular_momenta_;
  // The sum of the inertia tensors of the bodies with respect to their
  // individual centres of mass.  This is not the inertia tensor of the system
  // unless all bodies are at the same location, as their point masses also
  // contribute to the overall inertia.
  SymmetricBilinearForm<MomentOfInertia, InertialFrame> sum_of_inertia_tensors_;
};

}  // namespace internal_closed_system

using internal_closed_system::ClosedSystem;

}  // namespace physics
}  // namespace principia

#include "physics/closed_system_body.hpp"
