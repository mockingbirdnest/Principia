
#include "physics/euler_solver.hpp"

#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using quantities::Abs;
using quantities::Inverse;
using quantities::Sqrt;
using quantities::Square;
using quantities::SIUnit;
using quantities::si::Joule;
using quantities::si::Radian;

EulerSolver::EulerSolver(MomentOfInertia const& moment_of_inertia₁,
                         MomentOfInertia const& moment_of_inertia₂,
                         MomentOfInertia const& moment_of_inertia₃,
                         Energy const& kinetic_energy){
  // TODO(phl): What if they are not distinct?
  CHECK_LT(moment_of_inertia₁, moment_of_inertia₂);
  CHECK_LT(moment_of_inertia₂, moment_of_inertia₃);
  CHECK_LE(0 * Joule, kinetic_energy);

  AngularMomentum const G = SIUnit<AngularMomentum>();//TODO(phl):fixme.

  auto const& I₁ = moment_of_inertia₁;
  auto const& I₂ = moment_of_inertia₂;
  auto const& I₃ = moment_of_inertia₃;
  auto const T = kinetic_energy * (Radian * Radian);

  Square<AngularMomentum> const G² = G * G;
  auto const Δ₁ = G²  - 2.0 * T * I₁;
  auto const Δ₂ = G² - 2.0 * T * I₂;
  auto const Δ₃ = G² - 2.0 * T * I₃;

  auto const I₁₂ = Abs(I₁ - I₂);
  auto const I₁₃ = Abs(I₁ - I₃);
  auto const I₂₁ = Abs(I₂ - I₁);
  auto const I₂₃ = Abs(I₂ - I₃);
  auto const I₃₁ = Abs(I₃ - I₁);
  auto const I₃₂ = Abs(I₃ - I₂);

  B₁₃_ = Sqrt(I₁ * Δ₃ / I₁₃);
  B₃₁_ = Sqrt(I₃ * Δ₁ / I₃₁);

  // Note that Celledoni et al. give k, but we need mc = 1 - k^2.
  mc_ = 1.0 - Δ₁ * I₂₃ / (Δ₃ * I₂₁);
  λ₃_ = Sqrt(Δ₃ * I₁₂ / (I₁ * I₂ * I₃));

  //TODO(phl): Equality?
  if (2.0 * T * I₁ < G² && G² < 2.0 * T * I₂) {
    B₂₁_ = Sqrt(I₂ * Δ₁ / I₂₁);
    formula_ = Formula::i;
  } else if (2.0 * T * I₂ < G² && G² < 2.0 * T * I₃) {
    B₂₃_ = Sqrt(I₂ * Δ₃ / I₂₃);
    λ₁_ = Sqrt(Δ₁ * I₃₂ / (I₁ * I₂ * I₃));
    formula_ = Formula::ii;
  } else if (2.0 * T * I₂ == G²) {
    formula_ = Formula::iii;
  } else {
    LOG(FATAL) << "This should not happen";
  }
}

Bivector<AngularMomentum, EulerSolver::PrincipalAxesFrame>
EulerSolver::ComputeAngularMomentum(Instant const& t) {
  return Bivector<AngularMomentum, PrincipalAxesFrame>();
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
