
#include "physics/euler_solver.hpp"

#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using quantities::Abs;
using quantities::Sqrt;
using quantities::SIUnit;
using quantities::si::Joule;
using quantities::si::Radian;

EulerSolver::EulerSolver(MomentOfInertia const& moment_of_inertia₁,
                         MomentOfInertia const& moment_of_inertia₂,
                         MomentOfInertia const& moment_of_inertia₃,
                         Energy const& kinetic_energy)
    : moment_of_inertia₁_(moment_of_inertia₁),
      moment_of_inertia₂_(moment_of_inertia₂),
      moment_of_inertia₃_(moment_of_inertia₃),
      kinetic_energy_(kinetic_energy) {
  // TODO(phl): What if they are not distinct?
  CHECK_LT(moment_of_inertia₁_, moment_of_inertia₂_);
  CHECK_LT(moment_of_inertia₂_, moment_of_inertia₃_);
  CHECK_LE(0 * Joule, kinetic_energy_);

  AngularMomentum const G = SIUnit<AngularMomentum>();//TODO(phl):fixme.

  auto const& I₁ = moment_of_inertia₁;
  auto const& I₂ = moment_of_inertia₂;
  auto const& I₃ = moment_of_inertia₃;
  auto const T = Radian * Radian * kinetic_energy_ / (G * G);

  auto const Δ₁ = 1 - 2 * T * I₁;
  auto const Δ₂ = 1 - 2 * T * I₂;
  auto const Δ₃ = 1 - 2 * T * I₃;

  auto const I₁₂ = Abs(I₁ - I₂);
  auto const I₁₃ = Abs(I₁ - I₃);
  auto const I₂₁ = Abs(I₂ - I₁);
  auto const I₂₃ = Abs(I₂ - I₃);
  auto const I₃₁ = Abs(I₃ - I₁);
  auto const I₃₂ = Abs(I₃ - I₂);

  auto const B₁₃ = Sqrt(I₁ * Δ₃ / I₁₃);
  auto const B₂₁ = Sqrt(I₂ * Δ₁ / I₂₁);
  auto const B₂₃ = Sqrt(I₂ * Δ₃ / I₂₃);
  auto const B₃₁ = Sqrt(I₃ * Δ₁ / I₃₁);

  auto const k = Sqrt(Δ₁ * I₂₃ / (Δ₃ * I₂₁));//TODO(phl):We need mc!
  auto const λ₁ = Sqrt(Δ₁ * I₃₂ / (I₁ * I₂ * I₃));
  auto const λ₃ = Sqrt(Δ₃ * I₁₂ / (I₁ * I₂ * I₃));
}

Bivector<AngularMomentum, EulerSolver::PrincipalAxesFrame>
EulerSolver::ComputeAngularMomentum(Instant const& t) {
  return Bivector<AngularMomentum, PrincipalAxesFrame>();
}

}  // namespace internal_euler_solver
}  // namespace physics
}  // namespace principia
