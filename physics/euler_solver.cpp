
#include "physics/euler_solver.hpp"

#include <algorithm>

#include "numerics/elliptic_functions.hpp"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using geometry::Vector;
using numerics::EllipticF;
using numerics::JacobiSNCNDN;
using quantities::Abs;
using quantities::ArcTan;
using quantities::ArcTanh;
using quantities::Cosh;
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

EulerSolver::EulerSolver(
    R3Element<MomentOfInertia> const& moments_of_inertia,
    AngularMomentumBivector const& initial_angular_momentum,
    Instant const& initial_time)
    : initial_angular_momentum_(initial_angular_momentum),
      initial_time_(initial_time) {
  auto const& I₁ = moments_of_inertia.x;
  auto const& I₂ = moments_of_inertia.y;
  auto const& I₃ = moments_of_inertia.z;
  CHECK_LE(I₁, I₂);
  CHECK_LE(I₂, I₃);

  auto const& m = initial_angular_momentum.coordinates();

  auto const I₁₂ = I₁ - I₂;
  auto const I₁₃ = I₁ - I₃;
  auto const I₂₁ = -I₁₂;
  auto const I₂₃ = I₂ - I₃;
  auto const I₃₁ = -I₁₃;
  auto const I₃₂ = -I₂₃;

  // The formulæ for the Δs in Celledoni cannot be used directly because of
  // cancellations.
  auto const Δ₁ = m.y * m.y * I₂₁ / I₂ + m.z * m.z * I₃₁ / I₃;
  auto const Δ₂ = m.x * m.x * I₁₂ / I₁ + m.z * m.z * I₃₂ / I₃;
  auto const Δ₃ = m.x * m.x * I₃₁ / I₁ + m.y * m.y * I₃₂ / I₂;
  DCHECK_LE(Square<AngularMomentum>(), Δ₁);
  DCHECK_LE(Square<AngularMomentum>(), Δ₃);

  B₁₃_ = Sqrt(I₁ * Δ₃ / I₃₁);
  B₃₁_ = Sqrt(I₃ * Δ₁ / I₃₁);
  λ₃_ = Sqrt(Δ₃ * I₂₁ / (I₁ * I₂ * I₃));

  // Note that Celledoni et al. give k, but we need mc = 1 - k^2.  We write mc
  // in a way that reduces cancellations when k is close to 1.
  if (Δ₂ < Square<AngularMomentum>()) {
    B₂₁_ = Sqrt(I₂ * Δ₁ / I₂₁);
    mc_ = -Δ₂ * I₃₁ / (Δ₃ * I₂₁);
    ν_ = EllipticF(ArcTan(m.y / B₂₁_, m.z / B₃₁_), mc_) * Radian;
    if (m.x < AngularMomentum()) {
      λ₃_ = -λ₃_;
      B₁₃_ = -B₁₃_;
    }
    formula_ = Formula::i;
  } else if (Square<AngularMomentum>() < Δ₂) {
    B₂₃_ = Sqrt(I₂ * Δ₃ / I₃₂);
    mc_ = Δ₂ * I₃₁ / (Δ₁ * I₃₂);
    ν_ = EllipticF(ArcTan(m.y / B₂₃_, m.x / B₁₃_), mc_) * Radian;
    λ₁_ = Sqrt(Δ₁ * I₃₂ / (I₁ * I₂ * I₃));
    if (m.z < AngularMomentum()) {
      λ₁_ = -λ₁_;
      B₃₁_ = -B₃₁_;
    }
    formula_ = Formula::ii;
  } else {
    // Δ₂ == Square<AngularMomentum>()
    if (I₃₁ == MomentOfInertia()) {
      // The degenerate case of a sphere.  It would create NaNs.
      DCHECK_EQ(MomentOfInertia(), I₂₁);
      DCHECK_EQ(AngularFrequency(), λ₃_);
      formula_ = Formula::Sphere;
    } else {
      G_ =  initial_angular_momentum_.Norm();
      ν_ = -ArcTanh(m.y / G_);
      // NOTE(phl): The sign adjustments on this path are unclear.
      if (m.x < AngularMomentum()) {
        B₁₃_ = -B₁₃_;
      }
      if (m.z < AngularMomentum()) {
        B₃₁_ = -B₃₁_;
      }
      formula_ = Formula::iii;
    }
  }
}

EulerSolver::AngularMomentumBivector EulerSolver::AngularMomentumAt(
    Instant const& time) const {
  Time const Δt = time - initial_time_;
  switch (formula_) {
    case Formula::i: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN((λ₃_ * Δt - ν_) / Radian, mc_, sn, cn, dn);
      return AngularMomentumBivector({B₁₃_ * dn, -B₂₁_ * sn, B₃₁_ * cn});
    }
    case Formula::ii: {
      double sn;
      double cn;
      double dn;
      JacobiSNCNDN((λ₁_ * Δt - ν_) / Radian, mc_, sn, cn, dn);
      return AngularMomentumBivector({B₁₃_ * cn, -B₂₃_ * sn, B₃₁_ * dn});
    }
    case Formula::iii: {
      Angle const angle = λ₃_ * Δt - ν_;
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
      LOG(FATAL) << "Unexpected formula " << formula_;
  };
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
