#pragma once

#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {

using numerics::FixedVector;

namespace integrators {

// This class solves ordinary differential equations of following forms using a
// symplectic Runge-Kutta-Nyström method:
// (1).  (q, p)′ = X(q, p, t), with X = A(q, p) + B(q, p, t) and known
//       evolutions exp hA and exp hB, where [B, [B, [B, A]]] = 0;
// (2).  The above case, where (exp hB)(q, p, t) = (q, p) + h B(q, p, t),
//       (exp hA)(q, p) = (q, p) + h A(q, p), and A and B are known;
// (3).  q″ = f(q, t), with known f = -M⁻¹ ∇q V(q, t);
// Only (3) is implemented at this time.
// See the documentation for a proof that (3) is a special case of (2), and thus
// of (1), and for an explanation of the relation to Hamiltonian mechanics.

// Each step of size h is computed using the composition of evolutions
//   exp(aᵣ₋₁ h A) exp(bᵣ₋₁ h B) ... exp(a₀ h A)exp(b₀ h B);
// the integrator thus is a composition method.  If the appropriate coefficients
// vanish, the above can be reformulated as either
//   exp(aᵣ₋₁ h A) exp(bᵣ₋₁ h B) ... exp(b₁ h B) exp(a₀ h A)  or
//   exp(bᵣ₋₁ h B) exp(aᵣ₋₂ h A) ... exp(a₀ h A) exp(b₀ h B).
// The former is called type ABA, the latter type BAB, following the conventions
// used in Blanes, Casas and Ros (2001),
// New Families of Symplectic Runge-Kutta-Nyström Integration Methods,
// http://www.gicas.uji.es/Fernando/Proceedings/2000NAA.pdf.
// In the implementation, we call |stages_| the integer r above.  The number of
// |evaluations| is r-1 in the ABA and BAB cases, and r otherwise.
// See the documentation for an explanation of how types ABA and BAB reduce the
// number of evaluations required, especially in cases (2) and (3).

// As above, we follow the most widespread notation, calling the position
// weights a and the momentum weights b.  Note that our indices are 0-based.
// The following notations appear in the literature:
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

template<typename Position, int order_, bool time_reversible_, int evaluations_,
         CompositionMethod composition_>
class SymplecticRungeKuttaNyströmIntegrator
    : public FixedStepSizeIntegrator<
                 SpecialSecondOrderDifferentialEquation<Position>> {
  static int const stages_ = composition_ == kBA ? evaluations_
                                                 : evaluations_ + 1;
 public:
  SymplecticRungeKuttaNyströmIntegrator(FixedVector<double, stages_> const& a,
                                        FixedVector<double, stages_> const& b);

  void Solve(IntegrationProblem<ODE> const& problem,
             Time const& step) const override;

  static int const order = order_;
  static bool const time_reversible = time_reversible_;
  static int const evaluations = evaluations_;
  static CompositionMethod const composition = composition_;

 private:
  FixedVector<double, stages_> a_;
  FixedVector<double, stages_> b_;
  FixedVector<double, stages_> c_;
};

// This method minimizes the error constant.
// Coefficients from Robert I. McLachlan and Pau Atela (1992),
// The accuracy of symplectic integrators, table 2.
// http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      4 /*order*/,
                                      false /*time_reversible*/,
                                      4 /*evaluations*/,
                                      kBA> const&
McLachlanAtela1992Order4Optimal();
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      4 /*order*/,
                                      true /*time_reversible*/,
                                      4 /*evaluations*/,
                                      kABA> const& McLachlan1995SB3A4();
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      4 /*order*/,
                                      true /*time_reversible*/,
                                      5 /*evaluations*/,
                                      kABA> const& McLachlan1995SB3A5();
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      4 /*order*/,
                                      true /*time_reversible*/,
                                      6 /*evaluations*/,
                                      kBAB> const& BlanesMoan2002SRKN6B();
// This method minimizes the error constant.
// Coefficients from Robert I. McLachlan and Pau Atela (1992),
// The accuracy of symplectic integrators, table 2.
// http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      5 /*order*/,
                                      false /*time_reversible*/,
                                      6 /*evaluations*/,
                                      kBA> const&
McLachlanAtela1992Order5Optimal();
// Coefficients from Okunbor and Skeel (1994),
// Canonical Runge-Kutta-Nyström methods of orders 5 and 6,
// http://bionum.cs.purdue.edu/94OkSk.pdf.
// NOTE(egg): The coefficients were actually copied from McLachlan (1995), they
// seem to differ after a dozen significant figures or so.  Okunbor and Skeel
// remark "we did not use HYBRJ1 to improve the accuracy of method coefficients
// as we did in section 3.1".  We assume McLachlan's version is accurate.
// TODO(egg): Derive the coefficients with Mathematica.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      6 /*order*/,
                                      true /*time_reversible*/,
                                      7 /*evaluations*/,
                                      kABA> const&
OkunborSkeel1994Order6Method13();
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      6 /*order*/,
                                      true /*time_reversible*/,
                                      11 /*evaluations*/,
                                      kBAB> const& BlanesMoan2002SRKN11B();
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      6 /*order*/,
                                      true /*time_reversible*/,
                                      14 /*evaluations*/,
                                      kABA> const& BlanesMoan2002SRKN14A();

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_runge_kutta_nyström_integrator_body.hpp"
