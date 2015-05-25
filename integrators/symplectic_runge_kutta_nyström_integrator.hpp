#pragma once

#include "integrators/ordinary_differential_equations.hpp

namespace principia {

namespace integrators {

// This class solves ordinary differential equations of the forms
//   (q, p)′ = X, with X = A + B and known evolutions exp tA and exp tB,
//                where [B, [B, [B, A]]] = 0;       (1)
//   q″ = f(q, t), with known f = -M⁻¹ ∇q V(q);    (2)
// using a symplectic Runge-Kutta-Nyström method.
// Only (2) is implemented at this time.
// Note that (2) is a special case of (1), with
//   p = M q′, A = (M⁻¹ p, 0), B = (0, -∇q V(q)),
//   exp tA(q, p) = (q + t M⁻¹ p, p), exp tB(q, p) = (q, p - t ∇q V(q)).
// Also note that (2) corresponds to Hamilton's equations with
//   H(q, p) = ½ pᵀM⁻¹p + V(q).

// Each step of size h is computed using the composition of evolutions
//   exp(b₀ h B) exp(a₀ h A) ... exp(bᵣ h B) exp(aᵣ h A).
// If the appropriate coefficients vanish, this can be reformulated as either
//   exp(a₀ h A) exp(b₁ h B) ... exp(bᵣ h B) exp(aᵣ h A) or
//   exp(b₀ h B) exp(a₀ h A) ... exp(aᵣ₋₁ h A) exp(bᵣ h B).
// The former is called type ABA, the latter type BAB, following the conventions
// used in Blanes, Casas and Ros (2001), New families of symplectic
// Runge-Kutta-Nyström integration methods.

// Types ABA and BAB have the first-same-as-last property: the first and last
// applications of the evolution operators can be merged when output is not
// needed.

// For equations of type (2), methods of type BAB can make use of the FSAL
// property even for dense output, since d exp tB / dt is known, and need only
// be evaluated once for the two consecutive applications of exp tB.

enum CompositionFirstSameAsLast {
  kNone,
  kABA,
  kBAB,
};

template<typename Position, int order, int stages,
         CompositionFirstSameAsLast first_same_as_last>
class SymplecticRungeKuttaNyströmIntegrator
    : public FixedStepSizeIntegrator<
                 SpecialSecondOrderDifferentialEquation<Position>> {
  
  void Solve(IntegrationProblem<ODE> const& problem,
             Time const& step) const override;
};


}  // namespace integrators
}  // namespace principia
