
#include "physics/euler_solver.hpp"

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
    : initial_time_(initial_time) {
  auto const& I₁ = moments_of_inertia.x;
  auto const& I₂ = moments_of_inertia.y;
  auto const& I₃ = moments_of_inertia.z;

  // TODO(phl): What if they are not distinct?  We'll get infinities below and
  // that's probably fine.
  CHECK_LT(I₁, I₂);
  CHECK_LT(I₂, I₃);

  Square<AngularMomentum> const G² = initial_angular_momentum.Norm²();

  auto const& m = initial_angular_momentum.coordinates();
  auto const T = (m.x * m.x / I₁ + m.y * m.y / I₂ + m.z * m.z / I₃) * 0.5;

  auto const Δ₁ = Abs(G² - 2.0 * T * I₁);
  auto const Δ₂ = Abs(G² - 2.0 * T * I₂);
  auto const Δ₃ = Abs(G² - 2.0 * T * I₃);

  auto const I₁₂ = Abs(I₁ - I₂);
  auto const I₁₃ = Abs(I₁ - I₃);
  auto const I₂₁ = Abs(I₂ - I₁);
  auto const I₂₃ = Abs(I₂ - I₃);
  auto const I₃₁ = Abs(I₃ - I₁);
  auto const I₃₂ = Abs(I₃ - I₂);

  B₁₃_ = Sqrt(I₁ * Δ₃ / I₁₃);
  B₃₁_ = Sqrt(I₃ * Δ₁ / I₃₁);

  λ₃_ = Sqrt(Δ₃ * I₁₂ / (I₁ * I₂ * I₃));

  //TODO(phl): Equality?
  // Note that Celledoni et al. give k, but we need mc = 1 - k^2.
  if (2.0 * T * I₁ < G² && G² < 2.0 * T * I₂) {
    B₂₁_ = Sqrt(I₂ * Δ₁ / I₂₁);
    mc_ = 1.0 - Δ₁ * I₂₃ / (Δ₃ * I₂₁);
    ν_ = EllipticF(ArcTan(m.y / B₂₁_, m.z / B₃₁_), mc_) * Radian;
    if (m.x < AngularMomentum()) {
      λ₃_ = -λ₃_;
      B₁₃_ = -B₁₃_;
    }
    formula_ = Formula::i;
  } else if (2.0 * T * I₂ < G² && G² < 2.0 * T * I₃) {
    B₂₃_ = Sqrt(I₂ * Δ₃ / I₂₃);
    mc_ = 1.0 - Δ₃ * I₂₁ / (Δ₁ * I₂₃);
    ν_ = EllipticF(ArcTan(m.y / B₂₃_, m.x / B₁₃_), mc_) * Radian;
    λ₁_ = Sqrt(Δ₁ * I₃₂ / (I₁ * I₂ * I₃));
    if (m.z < AngularMomentum()) {
      λ₁_ = -λ₁_;
      B₃₁_ = -B₃₁_;
    }
    formula_ = Formula::ii;
  } else if (2.0 * T * I₂ == G²) {
    G_ =  Sqrt(G²);
    ν_= -ArcTanh(m.y / G_);
    // TODO(phl): The sign adjustments on this path are unclear.
    formula_ = Formula::iii;
  } else {
    LOG(FATAL) << "No formula for this case: G² = " << G²
               << ", 2.0 * T * I₁ = " << 2.0 * T * I₁
               << ", 2.0 * T * I₂ = " << 2.0 * T * I₂ << ", 2.0 * T * I₃ = "
               << 2.0 * T * I₃;
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
    default:
      LOG(FATAL) << "Unexpected formula";
  };
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
