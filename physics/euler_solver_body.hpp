
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
using geometry::DeduceSignPreservingOrientation;
using geometry::DefinesFrame;
using geometry::Normalize;
using geometry::Quaternion;
using geometry::Sign;
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
using quantities::IsFinite;
using quantities::Pow;
using quantities::Quotient;
using quantities::Sinh;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;
using quantities::Tanh;
using quantities::Time;
using quantities::Variation;
using quantities::si::Joule;
using quantities::si::Radian;

template<typename InertialFrame, typename PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>::EulerSolver(
    R3Element<MomentOfInertia> const& moments_of_inertia,
    Bivector<AngularMomentum, InertialFrame> const& initial_angular_momentum,
    AttitudeRotation const& initial_attitude,
    Instant const& initial_time)
    : moments_of_inertia_(moments_of_inertia),
      serialized_initial_angular_momentum_(initial_angular_momentum),
      initial_attitude_(initial_attitude),
      initial_time_(initial_time),
      G_(initial_angular_momentum.Norm()),
      ℛ_(Rotation<ℬʹ, InertialFrame>::Identity()),
      𝒮_(Signature<PrincipalAxesFrame,
                   PreferredPrincipalAxesFrame>::Identity()) {
  // Do not use initial_angular_momentum after this point.
  auto const initial_angular_momentum_in_principal_axes =
      initial_attitude.Inverse()(initial_angular_momentum);

  auto const& I₁ = moments_of_inertia_.x;
  auto const& I₂ = moments_of_inertia_.y;
  auto const& I₃ = moments_of_inertia_.z;
  CHECK_LE(I₁, I₂);
  CHECK_LE(I₂, I₃);

  // The usages of this variable prior to the computation of 𝒮_ must not depend
  // on the signs of its coordinates since we may flip it.
  auto m = initial_angular_momentum_in_principal_axes.coordinates();

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

  // These quantities are NaN in the spherical case, so they must be used with
  // care before we have checked for this case.
  auto const B₃₁² = I₃ * Δ₁ / I₃₁;
  auto const B₂₁² = I₂ * Δ₁ / I₂₁;
  auto const B₂₃² = I₂ * Δ₃ / I₂₃;
  auto const B₁₃² = I₁ * Δ₃ / I₁₃;
  B₁₃_ = Sqrt(B₁₃²);
  B₃₁_ = Sqrt(B₃₁²);

  // Determine the formula and region to use.
  if (Δ₂ < Square<AngularMomentum>()) {
    formula_ = Formula::i;
    region_ = Region::e₁;
  } else if (Square<AngularMomentum>() < Δ₂) {
    formula_ = Formula::ii;
    region_ = Region::e₃;
  } else {
    CHECK_EQ(Square<AngularMomentum>(), Δ₂);
    if (G_ == AngularMomentum()) {
      // No rotation.  Might as well handle it as a sphere.
      formula_ = Formula::Sphere;
      region_ = Region::Motionless;
    } else if (I₃₁ == MomentOfInertia()) {
      // The degenerate case of a sphere.  It would create NaNs.  Pick the
      // region that corresponds to the smallest coordinate.
      CHECK_EQ(MomentOfInertia(), I₂₁);
      CHECK_EQ(MomentOfInertia(), I₃₂);
      formula_ = Formula::Sphere;
      if (Abs(m.x) < Abs(m.z)) {
        region_ = Region::e₁;
      } else {
        region_ = Region::e₃;
      }
    } else {
      formula_ = Formula::iii;
      // Project along the smallest coordinate of x and z in absolute value.
      if (B₁₃_ < B₃₁_) {
        region_ = Region::e₁;
      } else {
        region_ = Region::e₃;
      }
    }
  }

  // Compute the rotation 𝒮_ that adjusts the signs of the coordinates of m in a
  // way that ensures that the quaternions are well-conditioned and that the σ's
  // disappear.
  Bivector<double, PreferredPrincipalAxesFrame> e₁({1, 0, 0});
  Bivector<double, PreferredPrincipalAxesFrame> e₂({0, 1, 0});
  Bivector<double, PreferredPrincipalAxesFrame> e₃({0, 0, 1});
  if (formula_ == Formula::iii) {
    𝒮_ = Signature<PrincipalAxesFrame, PreferredPrincipalAxesFrame>(
        Sign(m.x), DeduceSignPreservingOrientation{}, Sign(m.z));
  } else {
    switch (region_) {
      case Region::e₁: {
        𝒮_ = Signature<PrincipalAxesFrame, PreferredPrincipalAxesFrame>(
            Sign(m.x), Sign::Positive(), DeduceSignPreservingOrientation{});
        break;
      }
      case Region::e₃: {
        𝒮_ = Signature<PrincipalAxesFrame, PreferredPrincipalAxesFrame>(
            DeduceSignPreservingOrientation{}, Sign::Positive(), Sign(m.z));
        break;
      }
      case Region::Motionless: {
        𝒮_ = Signature<PrincipalAxesFrame,
                       PreferredPrincipalAxesFrame>::Identity();
        break;
      }
      default:
        LOG(FATAL) << "Unexpected region " << static_cast<int>(region_);
    }
  }

  // Now that 𝒮_ has been computed we can use it to adjust m and to compute ℛ_.
  initial_angular_momentum_ = 𝒮_(initial_angular_momentum_in_principal_axes);
  m = initial_angular_momentum_.coordinates();
  ℛ_ = [this, initial_attitude]() -> Rotation<ℬʹ, InertialFrame> {
    auto const 𝒴ₜ₀⁻¹ = Rotation<ℬʹ, ℬₜ>::Identity();
    auto const 𝒫ₜ₀⁻¹ = Compute𝒫ₜ(initial_angular_momentum_).Inverse();
    auto const 𝒮⁻¹ = 𝒮_.Inverse().template Forget<Rotation>();

    // This ℛ follows the assumptions in the third paragraph of section 2.3
    // of [CFSZ07], that is, the inertial frame is identified with the
    // (initial) principal axes frame.
    Rotation<ℬʹ, PrincipalAxesFrame> const ℛ = 𝒮⁻¹ * 𝒫ₜ₀⁻¹ * 𝒴ₜ₀⁻¹;

    // The multiplication by initial_attitude makes up for the loss of
    // generality due to the assumptions in the third paragraph of section
    // 2.3 of [CFSZ07].
    return initial_attitude * ℛ;
  }();

  switch (formula_) {
    case Formula::i: {
      CHECK_LE(Square<AngularMomentum>(), B₂₁²);
      B₂₁_ = Sqrt(B₂₁²);
      mc_ = std::min(Δ₂ * I₃₁ / (Δ₃ * I₂₁), 1.0);
      ν_ = EllipticF(ArcTan(m.y * B₃₁_, m.z * B₂₁_), mc_);
      auto const λ₃ = Sqrt(Δ₃ * I₁₂ / (I₁ * I₂ * I₃));
      λ_ = -λ₃;

      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(-ν_, mc_, sn, cn, dn);
      n_ = I₁ * I₃₂ / (I₃ * I₁₂);
      ψ_cn_multiplier_ = Sqrt(I₃ * I₂₁);
      ψ_sn_multiplier_ = Sqrt(I₂ * I₃₁);
      ψ_arctan_multiplier_ = -1;
      ψ_elliptic_pi_multiplier_ = G_ * I₁₃ / (λ_ * I₁ * I₃);
      ψ_offset_ = ψ_elliptic_pi_multiplier_ *
                      EllipticΠ(JacobiAmplitude(-ν_, mc_), n_, mc_) +
                  ψ_arctan_multiplier_ *
                      ArcTan(ψ_sn_multiplier_ * sn, ψ_cn_multiplier_ * cn);
      ψ_t_multiplier_ = G_ / I₁;

      break;
    }
    case Formula::ii: {
      CHECK_LE(Square<AngularMomentum>(), B₂₃²);
      B₂₃_ = Sqrt(B₂₃²);
      mc_ = std::min(Δ₂ * I₃₁ / (Δ₁ * I₃₂), 1.0);
      ν_ = EllipticF(ArcTan(m.y * B₁₃_, m.x * B₂₃_), mc_);
      auto const λ₁ = Sqrt(Δ₁ * I₃₂ / (I₁ * I₂ * I₃));
      λ_ = -λ₁;

      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(-ν_, mc_, sn, cn, dn);
      n_ = I₃ * I₂₁ / (I₁ * I₂₃);
      ψ_cn_multiplier_ = Sqrt(I₁ * I₃₂);
      ψ_sn_multiplier_ = Sqrt(I₂ * I₃₁);
      ψ_arctan_multiplier_ = 1;
      ψ_elliptic_pi_multiplier_ = G_ * I₃₁ / (λ_ * I₁ * I₃);
      ψ_offset_ = ψ_elliptic_pi_multiplier_ *
                      EllipticΠ(JacobiAmplitude(-ν_, mc_), n_, mc_) +
                  ψ_arctan_multiplier_ *
                      ArcTan(ψ_sn_multiplier_ * sn, ψ_cn_multiplier_ * cn);
      ψ_t_multiplier_ = G_ / I₃;

      break;
    }
    case Formula::iii: {
      ν_ = -ArcTanh(m.y / G_);
      auto const λ₂ = Sqrt(-Δ₁ * Δ₃ / (I₁ * I₃)) / G_;
      λ_ = λ₂;

      switch (region_) {
        case Region::e₁: {
          ψ_arctan_multiplier_ = -2;
          ψ_cosh_multiplier_ = B₃₁_;
          ψ_sinh_multiplier_ = B₁₃_ - G_;
          break;
        }
        case Region::e₃: {
          ψ_arctan_multiplier_ = 2;
          ψ_cosh_multiplier_ = B₁₃_;
          ψ_sinh_multiplier_ = B₃₁_ - G_;
          break;
        }
        case Region::Motionless:
        default:
          LOG(FATAL) << "Unexpected region " << static_cast<int>(region_);
      }
      ψ_offset_ = ArcTan(ψ_sinh_multiplier_ * Tanh(-0.5 * ν_),
                         ψ_cosh_multiplier_);
      ψ_t_multiplier_ = G_ / I₂;

      break;
    }
    case Formula::Sphere: {
      ψ_t_multiplier_ = G_ / I₂;
      break;
    }
  }
}

template<typename InertialFrame, typename PrincipalAxesFrame>
Bivector<AngularMomentum, PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>::AngularMomentumAt(
    Instant const& time) const {
  Time const Δt = time - initial_time_;
  PreferredAngularMomentumBivector m;
  switch (formula_) {
    case Formula::i: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      m = PreferredAngularMomentumBivector({B₁₃_ * dn, -B₂₁_ * sn, B₃₁_ * cn});
      break;
    }
    case Formula::ii: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      m = PreferredAngularMomentumBivector({B₁₃_ * cn, -B₂₃_ * sn, B₃₁_ * dn});
      break;
    }
    case Formula::iii: {
      Angle const angle = λ_ * Δt - ν_;
      double const sech = 1.0 / Cosh(angle);
      m = PreferredAngularMomentumBivector(
          {B₁₃_ * sech, G_ * Tanh(angle), B₃₁_ * sech});
      break;
    }
    case Formula::Sphere: {
      m = initial_angular_momentum_;
      break;
    }
    default:
      LOG(FATAL) << "Unexpected formula " << static_cast<int>(formula_);
  };
  return 𝒮_.Inverse()(m);
}

template<typename InertialFrame, typename PrincipalAxesFrame>
AngularVelocity<PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>::AngularVelocityFor(
    Bivector<AngularMomentum, PrincipalAxesFrame> const& angular_momentum)
    const {
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
    Bivector<AngularMomentum, PrincipalAxesFrame> const& angular_momentum,
    Instant const& time) const {
  Rotation<PreferredPrincipalAxesFrame, ℬₜ> const 𝒫ₜ =
      Compute𝒫ₜ(𝒮_(angular_momentum));

  Time const Δt = time - initial_time_;
  Angle ψ = ψ_t_multiplier_ * Δt;
  switch (formula_) {
    case Formula::i: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      Angle const φ = JacobiAmplitude(λ_ * Δt - ν_, mc_);
      ψ += ψ_elliptic_pi_multiplier_ * EllipticΠ(φ, n_, mc_) +
           ψ_arctan_multiplier_ *
               ArcTan(ψ_sn_multiplier_ * sn, ψ_cn_multiplier_ * cn) -
           ψ_offset_;
      break;
    }
    case Formula::ii: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN(λ_ * Δt - ν_, mc_, sn, cn, dn);
      Angle const φ = JacobiAmplitude(λ_ * Δt - ν_, mc_);
      ψ += ψ_elliptic_pi_multiplier_ * EllipticΠ(φ, n_, mc_) +
           ψ_arctan_multiplier_ *
               ArcTan(ψ_sn_multiplier_ * sn, ψ_cn_multiplier_ * cn) -
           ψ_offset_;
      break;
    }
    case Formula::iii: {
      ψ += ψ_arctan_multiplier_ *
           (ArcTan(ψ_sinh_multiplier_ * Tanh(0.5 * (λ_ * Δt - ν_)),
                   ψ_cosh_multiplier_) -
            ψ_offset_);
      break;
    }
    case Formula::Sphere: {
      break;
    }
    default:
      LOG(FATAL) << "Unexpected formula " << static_cast<int>(formula_);
  };

  switch (region_) {
    case Region::e₁: {
      Bivector<double, ℬʹ> const e₁({1, 0, 0});
      Rotation<ℬₜ, ℬʹ> const 𝒴ₜ(ψ, e₁, DefinesFrame<ℬₜ>{});
      return ℛ_ * 𝒴ₜ * 𝒫ₜ * 𝒮_.template Forget<Rotation>();
    }
    case Region::e₃: {
      Bivector<double, ℬʹ> const e₃({0, 0, 1});
      Rotation<ℬₜ, ℬʹ> const 𝒴ₜ(ψ, e₃, DefinesFrame<ℬₜ>{});
      return ℛ_ * 𝒴ₜ * 𝒫ₜ * 𝒮_.template Forget<Rotation>();
    }
    case Region::Motionless: {
      Bivector<double, ℬʹ> const unused({0, 1, 0});
      Rotation<ℬₜ, ℬʹ> const 𝒴ₜ(ψ, unused, DefinesFrame<ℬₜ>{});
      return ℛ_ * 𝒴ₜ * 𝒫ₜ * 𝒮_.template Forget<Rotation>();
    }
    default:
      LOG(FATAL) << "Unexpected region " << static_cast<int>(region_);
  }
}

template<typename InertialFrame, typename PrincipalAxesFrame>
R3Element<MomentOfInertia> const&
EulerSolver<InertialFrame, PrincipalAxesFrame>::moments_of_inertia() const {
  return moments_of_inertia_;
}

template<typename InertialFrame, typename PrincipalAxesFrame>
void EulerSolver<InertialFrame, PrincipalAxesFrame>::WriteToMessage(
    not_null<serialization::EulerSolver*> const message) const {
  moments_of_inertia_.WriteToMessage(message->mutable_moments_of_inertia());
  serialized_initial_angular_momentum_.WriteToMessage(
      message->mutable_initial_angular_momentum());
  initial_attitude_.WriteToMessage(message->mutable_initial_attitude());
  initial_time_.WriteToMessage(message->mutable_initial_time());
}

template<typename InertialFrame, typename PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>
EulerSolver<InertialFrame, PrincipalAxesFrame>::ReadFromMessage(
    serialization::EulerSolver const& message) {
  return EulerSolver(
      R3Element<MomentOfInertia>::ReadFromMessage(message.moments_of_inertia()),
      Bivector<AngularMomentum, InertialFrame>::ReadFromMessage(
          message.initial_angular_momentum()),
      AttitudeRotation::ReadFromMessage(message.initial_attitude()),
      Instant::ReadFromMessage(message.initial_time()));
}

template<typename InertialFrame, typename PrincipalAxesFrame>
Rotation<typename EulerSolver<InertialFrame,
                              PrincipalAxesFrame>::PreferredPrincipalAxesFrame,
         typename EulerSolver<InertialFrame, PrincipalAxesFrame>::ℬₜ>
EulerSolver<InertialFrame, PrincipalAxesFrame>::Compute𝒫ₜ(
    PreferredAngularMomentumBivector const& angular_momentum) const {
  auto const& m = angular_momentum;
  auto m_coordinates = m.coordinates();

  Quaternion pₜ;
  switch (region_) {
    case Region::e₁: {
      double const real_part = Sqrt(0.5 * (1 + m_coordinates.x / G_));
      AngularMomentum const denominator = 2 * G_ * real_part;
      pₜ = Quaternion(real_part,
                      {0,
                        m_coordinates.z / denominator,
                        -m_coordinates.y / denominator});
      break;
    }
    case Region::e₃: {
      double const real_part = Sqrt(0.5 * (1 + m_coordinates.z / G_));
      AngularMomentum const denominator = 2 * G_ * real_part;
      pₜ = Quaternion(real_part,
                      {m_coordinates.y / denominator,
                        -m_coordinates.x / denominator,
                        0});
      break;
    }
    case Region::Motionless: {
      pₜ = Quaternion(1);
      break;
    }
    default:
      LOG(FATAL) << "Unexpected region " << static_cast<int>(region_);
  }

  Rotation<PreferredPrincipalAxesFrame, ℬₜ> const 𝒫ₜ(pₜ);

  return 𝒫ₜ;
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
