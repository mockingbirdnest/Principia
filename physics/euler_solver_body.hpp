
#pragma once

#include "physics/euler_solver.hpp"

#include <algorithm>

#include "geometry/grassmann.hpp"
#include "geometry/quaternion.hpp"
#include "numerics/elliptic_functions.hpp"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using geometry::Commutator;
using geometry::DefinesFrame;
using geometry::Normalize;
using geometry::Quaternion;
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
EulerSolver<InertialFrame, PrincipalAxesFrame>::EulerSolver(
    R3Element<MomentOfInertia> const& moments_of_inertia,
    AngularMomentumBivector const& initial_angular_momentum,
    AttitudeRotation const& initial_attitude,
    Instant const& initial_time)
    : moments_of_inertia_(moments_of_inertia),
      initial_angular_momentum_(initial_angular_momentum),
      initial_time_(initial_time),
      ℛ_([this, initial_attitude]() -> Rotation<ℬʹ, InertialFrame> {
        auto const 𝒴ₜ₀⁻¹ = Rotation<ℬʹ, ℬₜ>::Identity();
        auto const 𝒫ₜ₀⁻¹ = Compute𝒫ₜ(initial_angular_momentum_).Inverse();

        // This ℛ follows the assumptions in the third paragraph of section 2.3
        // of [CFSZ07], that is, the inertial frame is identified with the
        // (initial) principal axes frame.
        Rotation<ℬʹ, PrincipalAxesFrame> const ℛ = 𝒫ₜ₀⁻¹ * 𝒴ₜ₀⁻¹;

        // The multiplication by initial_attitude makes up for the loss of
        // generality due to the assumptions in the third paragraph of section
        // 2.3 of [CFSZ07].
        return initial_attitude * ℛ;
      }()) {
  auto const& I₁ = moments_of_inertia_.x;
  auto const& I₂ = moments_of_inertia_.y;
  auto const& I₃ = moments_of_inertia_.z;
  CHECK_LE(I₁, I₂);
  CHECK_LE(I₂, I₃);

  auto const& m = initial_angular_momentum.coordinates();

  // These computations are such that if, say I₁ == I₂, I₂₁ is +0.0 and I₁₂ is
  // -0.0.
  auto const I₃₂ = I₃ - I₂;
  auto const I₃₁ = I₃ - I₁;
  auto const I₂₁ = I₂ - I₁;
  auto const I₂₃ = -I₃₂;
  auto const I₁₃ = -I₃₁;
  auto const I₁₂ = -I₂₁;

  // The formulæ for the Δs in [CFSZ07] cannot be used directly because of
  // cancellations.
  auto const Δ₁ = m.y * m.y * I₂₁ / I₂ + m.z * m.z * I₃₁ / I₃;
  auto const Δ₂ = m.z * m.z * I₃₂ / I₃ + m.x * m.x * I₁₂ / I₁;
  auto const Δ₃ = m.x * m.x * I₁₃ / I₁ + m.y * m.y * I₂₃ / I₂;
  CHECK_LE(Square<AngularMomentum>(), Δ₁);
  CHECK_LE(Δ₃, Square<AngularMomentum>());

  // These quantities are NaN in the spherical case, so they be used with care
  // before we have checked for this case.  We don't use B₁₂² and B₃₂² as they
  // might be negative with our definitions.
  auto const B₃₁² = I₃ * Δ₁ / I₃₁;
  auto const B₂₁² = I₂ * Δ₁ / I₂₁;
  auto const B₂₃² = I₂ * Δ₃ / I₂₃;
  auto const B₁₃² = I₁ * Δ₃ / I₁₃;
  B₁₃_ = Sqrt(B₁₃²);
  B₃₁_ = Sqrt(B₃₁²);

  auto const G² =  initial_angular_momentum_.Norm²();
  G_ =  Sqrt(G²);
  auto const two_T = m.x * m.x / I₁ + m.y * m.y / I₂ + m.z * m.z / I₃;
  ψ_t_multiplier_ = G_ / I₃;

  // Note that [CFSZ07] gives k, but we need mc = 1 - k^2.  We write mc in a way
  // that reduces cancellations when k is close to 1.
  if (Δ₂ < Square<AngularMomentum>()) {
    CHECK_LE(Square<AngularMomentum>(), B₂₃²);
    CHECK_LE(Square<AngularMomentum>(), B₂₁²);
    B₂₁_ = Sqrt(B₂₁²);
    mc_ = std::min(Δ₂ * I₃₁ / (Δ₃ * I₂₁), 1.0);
    ν_ = EllipticF(ArcTan(m.y * B₃₁_, m.z * B₂₁_), mc_);
    auto const λ₃ = Sqrt(Δ₃ * I₁₂ / (I₁ * I₂ * I₃));
    // TODO(phl): These tests on the signs of coordinates should probably handle
    // -0.0 correctly.
    if (m.x < AngularMomentum()) {
      σB₁₃_ = -B₁₃_;
      λ_ = λ₃;
    } else {
      σB₁₃_ = B₁₃_;
      λ_ = -λ₃;
    }

    double sn;
    double cn;
    double dn;
    JacobiSNCNDN(-ν_, mc_, sn, cn, dn);
    n_ = -B₃₁² / B₁₃²;
    ψ_arctan_x_multiplier_ = B₁₃_;
    ψ_arctan_y_multiplier_ = B₂₁_;
    ψ_arctan_multiplier_ = -B₃₁_ * ψ_arctan_x_multiplier_ /
                           (ψ_arctan_y_multiplier_ * G_);
    ψ_offset_ = EllipticΠ(JacobiAmplitude(-ν_, mc_), n_, mc_) +
                ψ_arctan_multiplier_ * ArcTan(ψ_arctan_y_multiplier_ * sn,
                                              ψ_arctan_x_multiplier_ * dn);
    ψ_integral_multiplier_ = G_ * I₃₁ / (λ_ * I₁ * I₃);

    formula_ = Formula::i;
  } else if (Square<AngularMomentum>() < Δ₂) {
    CHECK_LE(Square<AngularMomentum>(), B₂₃²);
    CHECK_LE(Square<AngularMomentum>(), B₂₁²);
    B₂₃_ = Sqrt(B₂₃²);
    mc_ = std::min(Δ₂ * I₃₁ / (Δ₁ * I₃₂), 1.0);
    ν_ = EllipticF(ArcTan(m.y * B₁₃_, m.x * B₂₃_), mc_);
    auto const λ₁ = Sqrt(Δ₁ * I₃₂ / (I₁ * I₂ * I₃));
    if (m.z < AngularMomentum()) {
      σB₃₁_ = -B₃₁_;
      λ_ = λ₁;
    } else {
      σB₃₁_ = B₃₁_;
      λ_ = -λ₁;
    }

    double sn;
    double cn;
    double dn;
    JacobiSNCNDN(-ν_, mc_, sn, cn, dn);
    n_ = I₂₁ * I₃ / (I₂₃ * I₁);
    ψ_arctan_x_multiplier_ = B₁₃_;
    ψ_arctan_y_multiplier_ = B₂₃_;
    ψ_arctan_multiplier_ = σB₃₁_ * ψ_arctan_x_multiplier_ /
                           (ψ_arctan_y_multiplier_ * G_);
    ψ_offset_ = EllipticΠ(JacobiAmplitude(-ν_, mc_), n_, mc_) +
                ψ_arctan_multiplier_ * ArcTan(ψ_arctan_y_multiplier_ * sn,
                                              ψ_arctan_x_multiplier_ * cn);
    ψ_integral_multiplier_ = G_ * I₃₁ / (λ_ * I₁ * I₃);

    formula_ = Formula::ii;
  } else {
    CHECK_EQ(Square<AngularMomentum>(), Δ₂);
    if (G_ == AngularMomentum()) {
      // No rotation.  Might as well handle it as a sphere.
      ψ_t_multiplier_ = AngularFrequency();
      formula_ = Formula::Sphere;
    } else if (I₃₁ == MomentOfInertia()) {
      // The degenerate case of a sphere.  It would create NaNs.
      CHECK_EQ(MomentOfInertia(), I₂₁);
      CHECK_EQ(MomentOfInertia(), I₃₂);
      formula_ = Formula::Sphere;
    } else {
      ν_ = -ArcTanh(m.y / G_);
      auto const λ₂ = Sqrt(-Δ₁ * Δ₃ / (I₁ * I₃)) / G_;
      λ_ = λ₂;
      if (m.x < AngularMomentum()) {
        σʹB₁₃_ = -B₁₃_;
        λ_ = -λ_;
      } else {
        σʹB₁₃_ = B₁₃_;
      }
      if (m.z < AngularMomentum()) {
        σʺB₃₁_ = -B₃₁_;
        λ_ = -λ_;
      } else {
        σʺB₃₁_ = B₃₁_;
      }
      // Δ₂ shows up in the multiplier, and λ_ is finite in the non-spherical
      // case, so things simplify tremendously.
      ψ_Π_offset_ = Angle();
      ψ_integral_multiplier_ = 0;
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
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      return AngularMomentumBivector({σB₁₃_ * dn, -B₂₁_ * sn, B₃₁_ * cn});
    }
    case Formula::ii: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      return AngularMomentumBivector({B₁₃_ * cn, -B₂₃_ * sn, σB₃₁_ * dn});
    }
    case Formula::iii: {
      Angle const angle = λ_ * Δt - ν_;
      double const sech = 1.0 / Cosh(angle);
      return AngularMomentumBivector(
          {σʹB₁₃_ * sech, G_ * Tanh(angle), σʺB₃₁_ * sech});
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
AngularVelocity<PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>::AngularVelocityFor(
    AngularMomentumBivector const& angular_momentum) const {
  auto const& m = angular_momentum;
  auto const& m_coordinates = m.coordinates();

  auto const& I₁ = moments_of_inertia_.x;
  auto const& I₂ = moments_of_inertia_.y;
  auto const& I₃ = moments_of_inertia_.z;
  Bivector<Quotient<AngularMomentum, MomentOfInertia>, PrincipalAxesFrame> const
      ω({m_coordinates.x / I₁, m_coordinates.y / I₂, m_coordinates.z / I₃});

  return ω;
}

template<typename InertialFrame, typename PrincipalAxesFrame>
typename EulerSolver<InertialFrame, PrincipalAxesFrame>::AttitudeRotation
EulerSolver<InertialFrame, PrincipalAxesFrame>::AttitudeAt(
    AngularMomentumBivector const& angular_momentum,
    Instant const& time) const {
  Rotation<PrincipalAxesFrame, ℬₜ> const 𝒫ₜ = Compute𝒫ₜ(angular_momentum);

  Time const Δt = time - initial_time_;
  Angle ψ = ψ_t_multiplier_ * Δt;
  switch (formula_) {
    case Formula::i: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      Angle const φ = JacobiAmplitude(λ_ * Δt - ν_, mc_);
      ψ += ψ_integral_multiplier_ *
           (EllipticΠ(φ, n_, mc_) +
            ψ_arctan_multiplier_ * ArcTan(ψ_arctan_y_multiplier_ * sn,
                                          ψ_arctan_x_multiplier_ * dn) -
            ψ_offset_);
      break;
    }
    case Formula::ii: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      Angle const φ = JacobiAmplitude(λ_ * Δt - ν_, mc_);
      ψ += ψ_integral_multiplier_ *
           (EllipticΠ(φ, n_, mc_) +
            ψ_arctan_multiplier_ * ArcTan(ψ_arctan_y_multiplier_ * sn,
                                          ψ_arctan_x_multiplier_ * cn) -
            ψ_offset_);
      break;
    }
    case Formula::iii: {
      break;
    }
    case Formula::Sphere: {
      break;
    }
    default:
      LOG(FATAL) << "Unexpected formula " << static_cast<int>(formula_);
  };
  Bivector<double, ℬʹ> const e₃({0, 0, 1});
  Rotation<ℬₜ, ℬʹ> const 𝒴ₜ(ψ, e₃, DefinesFrame<ℬₜ>{});

  return ℛ_ * 𝒴ₜ * 𝒫ₜ;
}

template<typename InertialFrame, typename PrincipalAxesFrame>
Rotation<PrincipalAxesFrame,
         typename EulerSolver<InertialFrame, PrincipalAxesFrame>::ℬₜ>
EulerSolver<InertialFrame, PrincipalAxesFrame>::Compute𝒫ₜ(
    AngularMomentumBivector const& angular_momentum) const {
  auto const& m = angular_momentum;
  auto m_coordinates = m.coordinates();

  // The first time through this function (at construction), determine if we'll
  // flip m.z and m.x to avoid the singularity m.z == -G.  After that, stick to
  // the decision made at construction.
  if (!must_flip_m_.has_value()) {
    must_flip_m_ = m_coordinates.z < AngularMomentum();
  }
  if (must_flip_m_.value()) {
    m_coordinates.x *= -1;
    m_coordinates.z *= -1;
  }

  double const real_part = Sqrt(0.5 * (1 + m_coordinates.z / G_));
  double const denominator = 2 * G_ * real_part;
  Quaternion const pₜ(real_part,
                      m_coordinates.y / denominator,
                      -m_coordinates.x / denominator,
                      0);

  Rotation<PrincipalAxesFrame, ℬₜ> const 𝒫ₜ(pₜ);

  return 𝒫ₜ;
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
