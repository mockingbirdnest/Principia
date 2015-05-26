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

// Note that (3) is a special case of (2), and thus of (1), with
//   p = M q′, A = (M⁻¹ p, 0), B = (0, -∇q V(q, t)),
//   (exp hA)(q, p) = (q + h M⁻¹ p, p), (exp hB)(q, p) = (q, p - h ∇q V(q)).
// It corresponds to Hamilton's equations with quadratic kinetic energy,
//   H(q, p, t) = ½ pᵀM⁻¹p + V(q, t).    (4)
// See below for a proof that these A and B do indeed satisfy
// [B, [B, [B, A]]] = 0.
//
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
// Types ABA and BAB have the first-same-as-last property: the last and first
// applications of the evolution operators in two consecutive steps can be
// merged when output is not needed, so for sparse outputs r-1 evolutions of B
// and r of A are required in the ABA case, and vice versa in the BAB case.
//
// When solving equations of type (2), the integrator makes full use of the FSAL
// property even for dense output: in the BAB case, the two consecutive
// applications of exp(bᵢ h B) require only one evaluation of B (and similarly
// in the ABA case).
//
// A remark on non-autonomy:
// Most treatments of these integrators write (1) as (q, p)′ = X, with
// X = A(q, p) + B(q, p), and (4) as H(q, p, t) = ½ pᵀM⁻¹p + V(q).
// It is however possible to incorporate time, by considering it as an
// additional variable:
//   (q, p, t)′ = X, with X = (A(q, p), 1) + (B(q, t, p), 0).
// Since B does not advance t, in cases of the form (2) we retain the defining
// property on the extended phase space,
//   (exp hB)(q, p, t) = (q, p, t) + (h B(q, p, t), 0).
//
// For equations of the form (3) it remains to show that Hamilton's equations
// with quadratic kinetic energy and a time-dependent potential satisfy
// [B, [B, [B, A]]] = 0.  Introducing t and its conjugate momentum ϖ to the
// phase space, (3) follows from Hamilton's equations with
//   H(q, p, t, ϖ) = ½ pᵀM⁻¹p + ϖ + V(q, t)
// since we then get t′ = 1.
// Writing Q = (q, t), P = (p, ϖ), and L(P) = ½ pᵀM⁻¹p + ϖ, the
// Hamiltonian becomes H(Q, P) = L(P) + V(Q), where L is a quadratic
// polynomial in P.
// Here A = {·, L}, and B = {·, V}, so that [B, [B, [B, A]]] = 0 is equivalent
// to {V, {V, {L, V}}} = 0, where {·, ·} is the Poisson bracket.  It is
// immediate every term in that expression will contain a third order partial
// derivative of L, and since L is quadratic all such derivatives vanish.  □
//
// See McLachlan and Quispel (2006), Geometric Integrators for ODEs, page 26,
// http://www.massey.ac.nz/~rmclachl/JPAReview.pdf for a detailed treatment
// of non-autonomous Hamiltonians using an extended phase space.
// See McLachlan (1993), Symplectic Integration of Wave Equations, page 8,
// http://www.massey.ac.nz/~rmclachl/wave.ps for a proof that
// {V, {V, {L, V}}} = 0 for arbitrary Poisson tensors and L quadratic in Q as
// well as P.

// Notations:
// As above, we the most widespread notation, calling the position weights aᵢ
// and the momentum weights bᵢ.  Note that our indices are 0-based.
// The following notations appear in the litterature:
//   (a, b) in most treatments of the subject;
//   (d, c) in Ruth, Yoshida, as well as Forest and Ruth;
//   (B, b) in Sofroniou and Spaletta;
//   (<unnamed>, B) in Okunbor and Skeel;
//   (<unnamed>, b) in Calvo and Sanz-Serna.
// Moreover, we follow the convention of Sofroniou and Spaletta (2002),
// Symplectic Methods for Separable Hamiltonian Systems, in calling cᵢ the
// nodes used for the time argument of the evolution of B, with
//   c₀ = 0, cᵢ = cᵢ₋₁ + aᵢ₋₁ for i > 0.
// The notation γᵢ is used for these nodes in Calvo and Sanz-Serna.
//
// Following common usage, we call the number of stages r - 1 in the ABA and BAB
// cases, and r otherwise.  Note that this is at odds with the traditions for
// nonsymplectic Runge-Kutta-Nyström methods, where the number of stages is
// considered to be the same regardless of the FSAL property.
//
// Note that the derivative weights b′ of general (nonsymplectic)
// Runge-Kutta-Nyström methods correspond to our momentum weights b; the
// position weights of the general case (usually denoted b) correspond to our
// bᵢ(1 - cᵢ), and the Runge-Kutta matrix is given by (usually denoted aᵤᵥ)
// is given by bᵥ (cᵤ - cᵥ).

enum CompositionKind {
  kNotFirstSameAsLast,
  kABA,
  kBAB,
};

template<typename Position, int order, int stages,
         CompositionKind first_same_as_last>
class SymplecticRungeKuttaNyströmIntegrator
    : public FixedStepSizeIntegrator<
                 SpecialSecondOrderDifferentialEquation<Position>> {
  
  void Solve(IntegrationProblem<ODE> const& problem,
             Time const& step) const override;
};


}  // namespace integrators
}  // namespace principia
