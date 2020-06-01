#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/quantities.hpp"
#include "physics/rigid_motion.hpp"

namespace principia {
namespace physics {
namespace internal_mechanical_system {

using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::InertiaTensor;
using geometry::Frame;
using geometry::OrthogonalMap;
using geometry::SymmetricBilinearForm;
using geometry::Vector;
using quantities::AngularMomentum;
using quantities::Mass;
using quantities::MomentOfInertia;

// Computes the instantaneous overall properties of a mechanical system.
// Effectively an extension of
// |BarycentreCalculator<DegreesOfFreedom<InertialFrame>, Mass>| that also
// handles rotational motion.
// |InertialFrame| is an inertial frame in which the motion of the member bodies
// is given; |SystemFrame| is a non-rotating frame with the same axes as
// |InertialFrame| and whose origin is the centre of mass of the system.
template<typename InertialFrame, typename SystemFrame>
class MechanicalSystem {
 public:
  static_assert(!InertialFrame::may_rotate);
  static_assert(!SystemFrame::may_rotate);

  MechanicalSystem() = default;

  template<typename BodyFrame>
  void AddRigidBody(RigidMotion<BodyFrame, InertialFrame> const& motion,
                    Mass const& mass,
                    InertiaTensor<BodyFrame> const& inertia_tensor);

  // The motion of a non-rotating frame whose origin is the centre of mass of
  // the system, and whose axes are those of |InertialFrame|.
  RigidMotion<SystemFrame, InertialFrame> LinearMotion() const;
  // The total mass.
  Mass const& mass() const;
  // The centre of mass.
  DegreesOfFreedom<InertialFrame> centre_of_mass() const;
  // The total angular momentum of the system, with respect to the origin of
  // |SystemFrame|, i.e., the centre of mass of the system.
  Bivector<AngularMomentum, SystemFrame> AngularMomentum() const;
  // The moments of inertia of the system, with respect to the origin of
  // |SystemFrame|, i.e., the centre of mass of the system.
  InertiaTensor<SystemFrame> InertiaTensor() const;

 private:
  BarycentreCalculator<DegreesOfFreedom<InertialFrame>, Mass> centre_of_mass_;
  std::vector<std::pair<DegreesOfFreedom<InertialFrame>, Mass>>
      body_linear_motions_;
  // The sum of the intrinsic angular momenta of the bodies, i.e., of the
  // angular momenta with respect to their individual centres of mass.  This is
  // not the total angular momentum, to which the linear motions of the bodies
  // with respect to the centre of mass of the system also contribute.
  Bivector<quantities::AngularMomentum, InertialFrame>
      sum_of_intrinsic_angular_momenta_;
  // The sum of the inertia tensors of the bodies with respect to their
  // individual centres of mass.  This is not the inertia tensor of the system
  // unless all bodies are at the same location, as their point masses also
  // contribute to the overall inertia.
  //
  // Note: It is more natural to manipulate the inertia tensor as a
  // |SymmetricBilinearForm<MomentOfInertia, Frame, Vector>| than a
  // |SymmetricBilinearForm<MomentOfInertia, Frame, Bivector>|.
  // The former has coordinates
  //   ⎛ ∑x²   ∑xy   ∑xy ⎞
  //   ⎜ ∑yx   ∑y²   ∑yz ⎟
  //   ⎝ ∑zx   ∑zy   ∑z² ⎠
  // whereas the latter (obtained from the former by |Anticommutator()|) has the
  // coordinates
  //   ⎛ ∑(y² + z²)   -∑xy      -∑xy    ⎞
  //   ⎜    -∑yx   ∑(x² + z²)   -∑yz    ⎟
  //   ⎝    -∑zx      -∑zy   ∑(x² + y²) ⎠
  // It is however common practice in physics to use the bilinear form on
  // bivectors, and to call its eigenvalues the principal moments of inertia.
  // This class thus exposes the form on bivectors in its API, and uses
  // |AnticommutatorInverse()| and |Anticommutator()| to convert between it and
  // the form on vectors used internally.
  SymmetricBilinearForm<MomentOfInertia, InertialFrame, Vector>
      sum_of_inertia_tensors_;
};

}  // namespace internal_mechanical_system

using internal_mechanical_system::MechanicalSystem;

}  // namespace physics
}  // namespace principia

#include "physics/mechanical_system_body.hpp"
