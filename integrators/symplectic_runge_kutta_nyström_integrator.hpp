#pragma once

#include "integrators/ordinary_differential_equations.hpp"

namespace principia {

namespace integrators {

// This class solves ordinary differential equations of following forms using a
// symplectic Runge-Kutta-Nyström method.
// (1).  (q, p)′ = X(q, p, t), with X = A(q, p) + B(q, p, t) and known
//       evolutions exp hA and exp hB, where [B, [B, [B, A]]] = 0;
// (2).  The above case, where (exp hB)(q, p, t) = (q, p) + h B(q, p, t),
//       (exp hA)(q, p) = (q, p) + h A(q, p), and A and B are known;
// (3).  q″ = f(q, t), with known f = -M⁻¹ ∇q V(q, t);
// Only (3) is implemented at this time.

// See the documentation for a proof that (3) is a special case of (2), and thus
// of (1), and for the relation to Hamiltonian mechanics.

// Each step of size h is computed using the composition of evolutions
//   exp(b₀ h B) exp(a₀ h A) ... exp(bᵣ h B) exp(aᵣ h A),
// so that the integrator is a composition method.  If the appropriate
// coefficients vanish, this can be reformulated as either
//   exp(a₀ h A) exp(b₁ h B) ... exp(bᵣ h B) exp(aᵣ h A) or
//   exp(b₀ h B) exp(a₀ h A) ... exp(aᵣ₋₁ h A) exp(bᵣ h B).
// The former is called type ABA, the latter type BAB, following the conventions
// used in Blanes, Casas and Ros (2001),
// New Families of Symplectic Runge-Kutta-Nyström Integration Methods,
// http://www.gicas.uji.es/Fernando/Proceedings/2000NAA.pdf.
// In the implementation, we call r the number of |stages_|.  The number of
// |evaluations| is r-1 in the ABA and BAB cases, and r otherwise.
// See the documentation for an explanation of how types ABA and BAB reduce the
// computational cost, especially in cases (2) and (3).

// As above, we follow the most widespread notation, calling the position
// weights a and the momentum weights b.  Note that our indices are 0-based.
// The following notations appear in the litterature:
//   (a, b) in most treatments of the subject;
//   (d, c) in Ruth, Yoshida, as well as Forest and Ruth;
//   (B, b) in Sofroniou and Spaletta;
//   (<unnamed>, B) in Okunbor and Skeel;
//   (<unnamed>, b) in Calvo and Sanz-Serna.
// Moreover, we follow the convention of Sofroniou and Spaletta (2002),
// Symplectic Methods for Separable Hamiltonian Systems, in calling c the
// nodes used for the time argument of the evolution of B, with
//   c₀ = 0, cᵢ = cᵢ₋₁ + aᵢ₋₁ for i > 0.
// The notation γ is used for these nodes in Calvo and Sanz-Serna.
// See the documentation for a description of the correspondence between
// these coefficients and those of a general Runge-Kutta-Nyström method.

enum CompositionMethod {
  kBA,   // Neither b₀ nor aᵣ vanishes.
  kABA,  // b₀ = 0.
  kBAB,  // aᵣ = 0.
};

template<typename Position, int order, int evaluations,
         CompositionMethod composition>
class SymplecticRungeKuttaNyströmIntegrator
    : public FixedStepSizeIntegrator<
                 SpecialSecondOrderDifferentialEquation<Position>> {
  static int const stages_ = composition == kBA ? evaluations : evaluations + 1;
 public:
  SymplecticRungeKuttaNyströmIntegrator(FixedVector<double, stages_> const& a,
                                        FixedVector<double, stages_> const& b,
                                        FixedVector<double, stages_> const& c);

  void Solve(IntegrationProblem<ODE> const& problem,
             Time const& step) const override;
 private:
  FixedVector<double, stages_> a_;
  FixedVector<double, stages_> b_;
  FixedVector<double, stages_> c_;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_runge_kutta_nyström_integrator_body.hpp"
