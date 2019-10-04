

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

using geometry::Commutator;
using geometry::Normalize;
using geometry::Vector;
using numerics::EllipticF;
using numerics::EllipticΠ;
using numerics::JacobiAmplitude;
using numerics::JacobiSNCNDN;
using quantities::Abs;
using quantities::ArcTan;
using quantities::ArcTanh;
using quantities::Cosh;
using quantities::Energy;
using quantities::Inverse;
using quantities::Pow;
using quantities::Quotient;
using quantities::Sinh;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;
using quantities::SIUnit;
using quantities::Tanh;
using quantities::Time;
using quantities::Variation;
using quantities::si::Joule;
using quantities::si::Radian;

template<typename InertialFrame, typename PrincipalAxesFrame>
Rotation<InertialFrame, PrincipalAxesFrame> Compute𝒫Rotation(
    R3Element<MomentOfInertia> const& moments_of_inertia,
    Bivector<AngularMomentum, PrincipalAxesFrame> const& angular_momentum) {
  auto const& m = angular_momentum;
  auto const& m_coordinates = m.coordinates();

  // Compute ṁ using the Euler equation.
  auto const& I₁ = moments_of_inertia.x;
  auto const& I₂ = moments_of_inertia.y;
  auto const& I₃ = moments_of_inertia.z;
  Bivector<Quotient<AngularMomentum, MomentOfInertia>, PrincipalAxesFrame> const
      ω({m_coordinates.x / I₁, m_coordinates.y / I₂, m_coordinates.z / I₃});
  Bivector<Variation<AngularMomentum>, PrincipalAxesFrame> const ṁ =
      Commutator(m, ω) / Radian;

  // Construct the orthonormal frame ℬₜ.
  auto const m_normalized = Normalize(m);
  auto const v = Normalize(ṁ);
  auto const w = Commutator(m_normalized, v);

  // 𝒫(m_normalized) = {0, 0, 1} , etc.
  Rotation<InertialFrame, PrincipalAxesFrame> const 𝒫(v, w, m_normalized);

  return 𝒫;
}

template<typename InertialFrame, typename PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>::EulerSolver(
    R3Element<MomentOfInertia> const& moments_of_inertia,
    AngularMomentumBivector const& initial_angular_momentum,
    AttitudeRotation const& initial_attitude,
    Instant const& initial_time)
    : moments_of_inertia_(moments_of_inertia),
      initial_angular_momentum_(initial_angular_momentum),
      initial_attitude_(initial_attitude),
      initial_time_(initial_time),
      𝒫ₜ₀_inverse_(Compute𝒫Rotation<InertialFrame, PrincipalAxesFrame>(
                       moments_of_inertia_,
                       initial_angular_momentum_)
                       .Inverse()) {
  auto const& I₁ = moments_of_inertia_.x;
  auto const& I₂ = moments_of_inertia_.y;
  auto const& I₃ = moments_of_inertia_.z;
  CHECK_LE(I₁, I₂);
  CHECK_LE(I₂, I₃);

  auto const& m = initial_angular_momentum.coordinates();

  auto const I₁₂ = I₁ - I₂;
  auto const I₂₃ = I₂ - I₃;
  auto const I₃₁ = I₃ - I₁;
  auto const I₂₁ = -I₁₂;
  auto const I₃₂ = -I₂₃;
  auto const I₁₃ = -I₃₁;

  // The formulæ for the Δs in Celledoni cannot be used directly because of
  // cancellations.
  auto const Δ₁ = m.y * m.y * I₂₁ / I₂ + m.z * m.z * I₃₁ / I₃;
  auto const Δ₂ = m.z * m.z * I₃₂ / I₃ + m.x * m.x * I₁₂ / I₁;
  auto const Δ₃ = m.x * m.x * I₁₃ / I₁ + m.y * m.y * I₂₃ / I₂;
  DCHECK_LE(Square<AngularMomentum>(), Δ₁);
  DCHECK_LE(Δ₃, Square<AngularMomentum>());

  auto const B₂₃² = I₂ * Δ₃ / I₂₃;
  auto const B₂₁² = I₂ * Δ₁ / I₂₁;
  DCHECK_LE(Square<AngularMomentum>(), B₂₃²);
  DCHECK_LE(Square<AngularMomentum>(), B₂₁²);

  B₁₃_ = Sqrt(I₁ * Δ₃ / I₁₃);
  B₃₁_ = Sqrt(I₃ * Δ₁ / I₃₁);

  auto const G² =  initial_angular_momentum_.Norm²();
  auto const two_T = m.x * m.x / I₁ + m.y * m.y / I₂ + m.z * m.z / I₃;
  ψ_t_multiplier_ = two_T / G_;

  // Note that Celledoni et al. give k, but we need mc = 1 - k^2.  We write mc
  // in a way that reduces cancellations when k is close to 1.
  if (Δ₂ < Square<AngularMomentum>()) {
    B₂₁_ = Sqrt(B₂₁²);
    mc_ = Δ₂ * I₃₁ / (Δ₃ * I₂₁);
    ν_ = EllipticF(ArcTan(m.y * B₃₁_, m.z * B₂₁_), mc_);
    λ₃_ = Sqrt(Δ₃ * I₁₂ / (I₁ * I₂ * I₃));
    if (m.x < AngularMomentum()) {
      B₁₃_ = -B₁₃_;
    } else {
      λ₃_ = -λ₃_;
    }
    n_ = G² / B₂₃²;
    ψ_Π_offset = EllipticΠ(-ν_, n_, mc_);
    ψ_Π_multiplier_ = Δ₂ / (λ₃_ * I₂ * G_);
    formula_ = Formula::i;
  } else if (Square<AngularMomentum>() < Δ₂) {
    B₂₃_ = Sqrt(B₂₃²);
    mc_ = Δ₂ * I₃₁ / (Δ₁ * I₃₂);
    ν_ = EllipticF(ArcTan(m.y * B₁₃_, m.x * B₂₃_), mc_);
    λ₁_ = Sqrt(Δ₁ * I₃₂ / (I₁ * I₂ * I₃));
    if (m.z < AngularMomentum()) {
      B₃₁_ = -B₃₁_;
    } else {
      λ₁_ = -λ₁_;
    }
    n_ = G² / B₂₁²;
    ψ_Π_offset = EllipticΠ(-ν_, n_, mc_);
    ψ_Π_multiplier_ = Δ₂ / (λ₁_ * I₂ * G_);
    formula_ = Formula::ii;
  } else {
    CHECK_EQ(Square<AngularMomentum>(), Δ₂);
    if (I₃₁ == MomentOfInertia()) {
      // The degenerate case of a sphere.  It would create NaNs.
      DCHECK_EQ(MomentOfInertia(), I₂₁);
      DCHECK_EQ(MomentOfInertia(), I₃₂);
      formula_ = Formula::Sphere;
    } else {
      G_ =  Sqrt(G²);
      ν_ = -ArcTanh(m.y / G_);
      λ₂_ = Sqrt(-Δ₁ * Δ₃ / (I₁ * I₃)) / G_;
      if (m.x < AngularMomentum()) {
        B₁₃_ = -B₁₃_;
        λ₂_ = -λ₂_;
      }
      if (m.z < AngularMomentum()) {
        B₃₁_ = -B₃₁_;
        λ₂_ = -λ₂_;
      }
      // Not quite an elliptic integral characteristic, but we'll stick to that
      // notation.
      n_ = G² * G² / (B₂₁² * B₂₃²);
      ψ_Π_offset = (-ν_ + n_ * std::log(n_ * Sinh(-ν_) - Cosh(-ν_)) * Radian);
      ψ_Π_multiplier_ = Δ₂ / (λ₂_ * I₂ * G_ * (1 - n_ * n_));
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
  Rotation<InertialFrame, PrincipalAxesFrame> const 𝒫ₜ =
      Compute𝒫Rotation<InertialFrame, PrincipalAxesFrame>(moments_of_inertia_,
                                                          angular_momentum);

  Time const Δt = time - initial_time_;
  Angle ψ;
  switch (formula_) {
    case Formula::i: {
      Angle const φ = JacobiAmplitude(λ₃_ * Δt - ν_, mc_);
      ψ = ψ_t_multiplier_ * Δt +
          ψ_Π_multiplier_ * (EllipticΠ(φ, n_, mc_) - ψ_Π_offset);
      break;
    }
    case Formula::ii: {
      Angle const φ = JacobiAmplitude(λ₁_ * Δt - ν_, mc_);
      ψ = ψ_t_multiplier_ * Δt +
          ψ_Π_multiplier_ * (EllipticΠ(φ, n_, mc_) - ψ_Π_offset);
      break;
    }
    case Formula::iii: {
      Angle const angle = λ₂_ * Δt - ν_;
      ψ = ψ_t_multiplier_ * Δt +
          ψ_Π_multiplier_ *
              (angle + n_ * std::log(n_ * Sinh(angle) - Cosh(angle)) * Radian -
               ψ_Π_offset);
      break;
    }
    case Formula::Sphere: {
    }
    default:
      LOG(FATAL) << "Unexpected formula " << static_cast<int>(formula_);
  };
  Bivector<double, PrincipalAxesFrame> const e₃({0, 0, 1});
  Rotation<PrincipalAxesFrame, PrincipalAxesFrame> const 𝒴ₜ(ψ, e₃);

  return 𝒫ₜ₀_inverse_ * 𝒴ₜ * 𝒫ₜ;
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
