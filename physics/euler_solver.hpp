
#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using geometry::Bivector;
using geometry::Instant;
using geometry::R3Element;
using geometry::Rotation;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::AngularMomentum;
using quantities::MomentOfInertia;
using quantities::NaN;

// A solver for Euler's rotation equations.  It follows Celledoni, Fassò,
// Säfström and Zanna, 2007, The exact computation of the free rigid body motion
// and its use in splitting method.
// NOTE(phl): There are a number of errors in the formulæ in Proposition 2.1, as
// can be seen by differentiation:
//   In case (i) λ should be defined as -σ λ₃.
//   In case (ii) λ should be defined as -σ λ₁.
//   In case (iii) the first coordinate should include a factor σ (not σʹ) and
//   λ should be defined as σ σʹ λ₂ (where λ₂ is the common value of λ₁ and λ₃).
template<typename InertialFrame, typename PrincipalAxesFrame>
class EulerSolver {
  static_assert(InertialFrame::is_inertial);
 public:
  using AngularMomentumBivector = Bivector<AngularMomentum, PrincipalAxesFrame>;
  using AttitudeRotation = Rotation<InertialFrame, PrincipalAxesFrame>;

  // Constructs a solver for a body with the given moments_of_inertia in its
  // principal axes frame.  The moments must be in increasing order.  At
  // initial_time the angular momentum is initial_angular_momentum and the
  // attitude initial_attitude.
  EulerSolver(R3Element<MomentOfInertia> const& moments_of_inertia,
              AngularMomentumBivector const& initial_angular_momentum,
              AttitudeRotation const& initial_attitude,
              Instant const& initial_time);

  // Computes the angular momentum at the given time.
  AngularMomentumBivector AngularMomentumAt(Instant const& time) const;

  // Computes the attitude at the given time, using the angular momentum
  // computed by the previous function.
  AttitudeRotation AttitudeAt(AngularMomentumBivector const& angular_momentum,
                              Instant const& time) const;

 private:
  // The formula to use, following Cellodoni et al., Section 2.2.  They don't
  // have a formula for the spherical case.
  enum class Formula {
    i,
    ii,
    iii,
    Sphere
  };

  // Construction parameters.
  R3Element<MomentOfInertia> const moments_of_inertia_;
  AngularMomentumBivector const initial_angular_momentum_;
  AttitudeRotation const initial_attitude_;
  Instant const initial_time_;

  MomentOfInertia const& I₁_ = moments_of_inertia_.x;
  MomentOfInertia const& I₂_ = moments_of_inertia_.y;
  MomentOfInertia const& I₃_ = moments_of_inertia_.z;

  // Amusingly, the formula to use is a constant of motion.
  Formula formula_;

  // Only the parameters needed for the selected formula are non-NaN after
  // construction.
  AngularMomentum B₁₃_ = NaN<AngularMomentum>();
  AngularMomentum B₃₁_ = NaN<AngularMomentum>();
  AngularMomentum B₂₁_ = NaN<AngularMomentum>();
  AngularMomentum B₂₃_ = NaN<AngularMomentum>();
  AngularMomentum G_ = NaN<AngularMomentum>();
  AngularFrequency λ₁_ = NaN<AngularFrequency>();
  AngularFrequency λ₂_ = NaN<AngularFrequency>();
  AngularFrequency λ₃_ = NaN<AngularFrequency>();
  double mc_ = NaN<double>();
  Angle ν_ = NaN<Angle>();
};

}  // namespace internal_euler_solver

using internal_euler_solver::EulerSolver;

}  // namespace physics
}  // namespace principia

#include "physics/euler_solver_body.hpp"
