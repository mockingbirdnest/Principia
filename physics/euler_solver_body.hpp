
#pragma once

#include "physics/euler_solver.hpp"

#include <algorithm>

#include "geometry/grassmann.hpp"
#include "numerics/elliptic_functions.hpp"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using geometry::Vector;
using numerics::EllipticF;
using numerics::EllipticΠ;
using numerics::JacobiAmplitude;
using numerics::JacobiSNCNDN;
using quantities::Abs;
using quantities::ArcTan;
using quantities::ArcTanh;
using quantities::Cosh;
using quantities::Pow;
using quantities::Tanh;
using quantities::Energy;
using quantities::Inverse;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;
using quantities::SIUnit;
using quantities::Time;
using quantities::si::Joule;
using quantities::si::Radian;

template<typename InertialFrame, typename PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>::EulerSolver(
    R3Element<MomentOfInertia> const& moments_of_inertia,
    AngularMomentumBivector const& initial_angular_momentum,
    AttitudeRotation const& initial_attitude,
    Instant const& initial_time)
    : moments_of_inertia_(moments_of_inertia),
      initial_angular_momentum_(initial_angular_momentum),
      initial_attitude_(initial_attitude),
      initial_time_(initial_time) {
  CHECK_LE(I₁_, I₂_);
  CHECK_LE(I₂_, I₃_);

  auto const& m = initial_angular_momentum.coordinates();

  auto const I₁₂ = I₁_ - I₂_;
  auto const I₁₃ = I₁_ - I₃_;
  auto const I₂₁ = -I₁₂;
  auto const I₂₃ = I₂_ - I₃_;
  auto const I₃₁ = -I₁₃;
  auto const I₃₂ = -I₂₃;

  // The formulæ for the Δs in Celledoni cannot be used directly because of
  // cancellations.
  auto const Δ₁ = m.y * m.y * I₂₁ / I₂_ + m.z * m.z * I₃₁ / I₃_;
  auto const Δ₂ = m.x * m.x * I₁₂ / I₁_ + m.z * m.z * I₃₂ / I₃_;
  auto const Δ₃ = m.x * m.x * I₃₁ / I₁_ + m.y * m.y * I₃₂ / I₂_;
  DCHECK_LE(Square<AngularMomentum>(), Δ₁);
  DCHECK_LE(Square<AngularMomentum>(), Δ₃);

  B₁₃_ = Sqrt(I₁_ * Δ₃ / I₃₁);
  B₃₁_ = Sqrt(I₃_ * Δ₁ / I₃₁);

  // Note that Celledoni et al. give k, but we need mc = 1 - k^2.  We write mc
  // in a way that reduces cancellations when k is close to 1.
  if (Δ₂ < Square<AngularMomentum>()) {
    B₂₁_ = Sqrt(I₂_ * Δ₁ / I₂₁);
    n_ = 1 / B₂₃_;  //?
    mc_ = -Δ₂ * I₃₁ / (Δ₃ * I₂₁);
    ν_ = EllipticF(ArcTan(m.y / B₂₁_, m.z / B₃₁_), mc_);
    λ₃_ = Sqrt(Δ₃ * I₂₁ / (I₁_ * I₂_ * I₃_));
    if (m.x < AngularMomentum()) {
      B₁₃_ = -B₁₃_;
    } else {
      λ₃_ = -λ₃_;
    }
    ψ_offset_ = EllipticΠ(-ν_, n_, mc_);
    ψ_multiplier_ = Δ₂ / (λ₃_ * I₂_);  //?
    formula_ = Formula::i;
  } else if (Square<AngularMomentum>() < Δ₂) {
    B₂₃_ = Sqrt(I₂_ * Δ₃ / I₃₂);
    mc_ = Δ₂ * I₃₁ / (Δ₁ * I₃₂);
    ν_ = EllipticF(ArcTan(m.y / B₂₃_, m.x / B₁₃_), mc_);
    λ₁_ = Sqrt(Δ₁ * I₃₂ / (I₁_ * I₂_ * I₃_));
    if (m.z < AngularMomentum()) {
      B₃₁_ = -B₃₁_;
    } else {
      λ₁_ = -λ₁_;
    }
    formula_ = Formula::ii;
  } else {
    CHECK_EQ(Square<AngularMomentum>(), Δ₂);
    if (I₃₁ == MomentOfInertia()) {
      // The degenerate case of a sphere.  It would create NaNs.
      DCHECK_EQ(MomentOfInertia(), I₂₁);
      DCHECK_EQ(MomentOfInertia(), I₃₂);
      formula_ = Formula::Sphere;
    } else {
      G_ =  initial_angular_momentum_.Norm();
      ν_ = -ArcTanh(m.y / G_);
      λ₂_ = Sqrt(Δ₁ * Δ₃ / (I₁_ * I₃_)) / G_;
      if (m.x < AngularMomentum()) {
        B₁₃_ = -B₁₃_;
        λ₂_ = -λ₂_;
      }
      if (m.z < AngularMomentum()) {
        B₃₁_ = -B₃₁_;
        λ₂_ = -λ₂_;
      }
      formula_ = Formula::iii;
    }
  }
}

template<typename InertialFrame, typename PrincipalAxesFrame>
typename EulerSolver<InertialFrame, PrincipalAxesFrame>::AngularMomentumBivector
EulerSolver<InertialFrame, PrincipalAxesFrame>::AngularMomentumAt(
    Instant const& time) const {
  Time const Δt = time - initial_time_;
  switch (formula_) {
    case Formula::i: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ₃_ * Δt - ν_, mc_, sn, cn, dn);
      return AngularMomentumBivector({B₁₃_ * dn, -B₂₁_ * sn, B₃₁_ * cn});
    }
    case Formula::ii: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ₁_ * Δt - ν_, mc_, sn, cn, dn);
      return AngularMomentumBivector({B₁₃_ * cn, -B₂₃_ * sn, B₃₁_ * dn});
    }
    case Formula::iii: {
      Angle const angle = λ₂_ * Δt - ν_;
      double const sech = 1.0 / Cosh(angle);
      return AngularMomentumBivector(
          {B₁₃_ * sech, G_ * Tanh(angle), B₃₁_ * sech});
    }
    case Formula::Sphere : {
      // NOTE(phl): It's unclear how the formulæ degenerate in this case, but
      // surely λ₃_ becomes 0, so the dependency in time disappears, so this is
      // my best guess.
      return initial_angular_momentum_;
    }
    default:
      LOG(FATAL) << "Unexpected formula " << static_cast<int>(formula_);
  };
}

template<typename InertialFrame, typename PrincipalAxesFrame>
typename EulerSolver<InertialFrame, PrincipalAxesFrame>::AttitudeRotation
EulerSolver<InertialFrame, PrincipalAxesFrame>::AttitudeAt(
    AngularMomentumBivector const& angular_momentum,
    Instant const& time) const {
  auto const& m = angular_momentum.coordinates();
  Time const Δt = time - initial_time_;
  switch (formula_) {
    case Formula::i: {
      // Note that the sign of λ has been integrated in λ₃_ at construction.
      Angle const φ = JacobiAmplitude(λ₃_ * Δt - ν_, mc_);
      Angle const ψ =
          2 * T_ * Δt + ψ_multiplier_ * (EllipticΠ(φ, n_, mc_) - ψ_offset_);
    }
    case Formula::ii: {
    }
    case Formula::iii: {
    }
    case Formula::Sphere : {
    }
    default:
      LOG(FATAL) << "Unexpected formula " << static_cast<int>(formula_);
  };
  return AttitudeRotation::Identity();
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
