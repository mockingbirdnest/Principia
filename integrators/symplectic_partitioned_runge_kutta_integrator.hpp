
#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_integrator.hpp"

namespace principia {

using base::not_null;

namespace integrators {

// The integrator defined in this class integrates Hamilton's equations for
// Hamiltonians of the form
//   H(q, p, t) = T(p) + V(q, t),
// where ∂T/∂pᵢ(q) and ∂V/∂qᵢ(q,t) are known.
// The asymmetric support for time dependence arises in the treatment of t as a
// variable in the extended phase space (q, t, p, ϖ) with the extended
// Hamiltonian  H(q, t, p, ϖ) = (T(p) + ϖ) + V(q, t), which can then be 
// integrated with methods for time-independent separable Hamiltonians
// H(Q, P) = T(P) + V(Q).
// For practical purposes, this reduces to evaluating the force computations
// with an additional parameter given by the weights (in the usual Runge-Kutta
// sense) corresponding to the positions.
// See McLachlan and Quispel (2006), Geometric Integrators for ODEs,
// http://www.massey.ac.nz/~rmclachl/JPAReview.pdf,
// or Sofroniou and Spaletta (2002),
// Symplectic Methods for Separable Hamiltonian Systems.

// We follow the convention of McLachlan & Atela, calling the position nodes
// aᵢ and the momentum nodes bᵢ.  The following notations appear in the
// litterature:
//   (a, b) in McLachlan & Atela, as well as Candy & Rozmus;
//   (d, c) in Ruth, Yoshida, as well as Forest & Ruth;
//   (B, b) in Sofroniou & Spaletta.
// Moreover, we follow the convention of Sofroniou & Spaletta in calling the
// weights used for the time argument of the force computation cᵢ, with
// c₁ = 0, cᵢ = cᵢ₋₁ + aᵢ₋₁ for i > 1.

enum VanishingCoefficients {
  kNone,
  kFirstBVanishes,
  kLastAVanishes,
};

template<typename Position, typename Momentum>
class SPRKIntegrator : public SymplecticIntegrator<Position, Momentum> {
 public:
  using Parameters = typename SymplecticIntegrator<Position,
                                                   Momentum>::Parameters;
  using Scheme = typename SymplecticIntegrator<Position, Momentum>::Scheme;
  using SystemState = typename SymplecticIntegrator<Position,
                                                    Momentum>::SystemState;

  SPRKIntegrator();
  ~SPRKIntegrator() override = default;

  // First-same-as-last (FSAL) methods with k stages require k - 1 force
  // and velocity evaluations for sparse output.  For dense output, methods
  // with synchronous momenta require k - 1 force evaluations and k velocity
  // evaluations, methods with synchronous positions require k force evaluations
  // and k - 1 velocity evaluations.

  // Second order, 2 stages, FSAL (synchronous momenta).
  Scheme const& Leapfrog() const;
  // Second order, 2 stages, FSAL (synchronous positions).
  Scheme const& PseudoLeapfrog() const;
  // Second order, 2 stages.  This method minimizes the error constant.
  // Coefficients from McLachlan and Atela (1992),
  // The accuracy of symplectic integrators, table 2.
  // http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
  Scheme const& McLachlanAtela1992Order2Optimal() const;
  // Third order, 3 stages.
  // Coefficients from Ruth (1983), A canonical integration technique,
  // https://accelconf.web.cern.ch/accelconf/p83/PDF/PAC1983_2669.PDF.
  Scheme const& Ruth1983() const;
  // Third order, 3 stages.  This method minimizes the error constant.
  // Coefficients from McLachlan and Atela (1992),
  // The accuracy of symplectic integrators, equation 3.3.
  // http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
  // NOTE(egg): the coefficients given in table 2 for this integrator are
  // incorrect (b1 and b2 are swapped).
  Scheme const& McLachlanAtela1992Order3Optimal() const;
  // Fourth order, 4 stages, FSAL (synchronous momenta).
  // Coefficients from Forest and Ruth (1990),
  // Fourth-order symplectic integration, equation 4.8.
  // http://zwe.web.cern.ch/zwe/CAS/biblio/ruth-forest.pdf.
  // This scheme was independently discovered by Candy and Rozmus (1991),
  // A Symplectic Integration Algorithm for Separable Hamiltonian Functions
  // (submitted earlier and published later than the Forest and Ruth paper).
  Scheme const& CandyRozmus1991ForestRuth1990SynchronousMomenta() const;
  // Fourth order, 4 stages, FSAL (synchronous positions).
  // Coefficients from Forest and Ruth (1990),
  // Fourth-order symplectic integration, equation 4.9.
  // http://zwe.web.cern.ch/zwe/CAS/biblio/ruth-forest.pdf.
  // This scheme was independently discovered by Candy and Rozmus (1991),
  // A Symplectic Integration Algorithm for Separable Hamiltonian Functions
  // (submitted earlier and published later than the Forest and Ruth paper).
  Scheme const& CandyRozmus1991ForestRuth1990SynchronousPositions() const;
  // Fourth order, 4 stages.  This method minimizes the error constant.
  // Coefficients from Robert I. McLachlan and Pau Atela (1992),
  // The accuracy of symplectic integrators, table 2.
  // http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
  Scheme const& McLachlanAtela1992Order4Optimal() const;
  // Fifth order, 6 stages.  This method minimizes the error constant  Ibidem.
  Scheme const& McLachlanAtela1992Order5Optimal() const;
  // Sixth order, 8 stages, FSAL (synchronous positions).
  // Coefficients from Yoshida (1990),
  // Construction of higher order symplectic integrators
  // http://sixtrack.web.cern.ch/SixTrack/doc/yoshida00.pdf.
  // NOTE(egg): The coefficients were derived from equations 5.4 through 5.17
  // rather than computed from the wᵢ given in tables 1 and 2.  The results were
  // then cross-checked agains those obtained from the tables.
  Scheme const& Yoshida1990Order6A() const;
  // Sixth order, 8 stages, FSAL (synchronous positions).  Ibidem.
  Scheme const& Yoshida1990Order6B() const;
  // Sixth order, 8 stages, FSAL (synchronous positions).  Ibidem.
  Scheme const& Yoshida1990Order6C() const;
  // Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
  Scheme const& Yoshida1990Order8A() const;
  // Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
  Scheme const& Yoshida1990Order8B() const;
  // Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
  Scheme const& Yoshida1990Order8C() const;
  // Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
  Scheme const& Yoshida1990Order8D() const;
  // Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
  Scheme const& Yoshida1990Order8E() const;

  void Initialize(Scheme const& coefficients) override;

  // The functors |compute_velocity| and |compute_force| compute
  // ∂T/∂pᵢ(p) and ∂V/∂qᵢ(q,t) respectively.
  template<typename AutonomousRightHandSideComputation,
           typename RightHandSideComputation>
  void Solve(RightHandSideComputation compute_force,
             AutonomousRightHandSideComputation compute_velocity,
             Parameters const& parameters,
             not_null<std::vector<SystemState>*> const solution) const;

 private:
  struct FirstSameAsLast {
    double first;
    double last;
  };

  template<VanishingCoefficients vanishing_coefficients,
           typename AutonomousRightHandSideComputation,
           typename RightHandSideComputation>
  void SolveOptimized(RightHandSideComputation compute_force,
                      AutonomousRightHandSideComputation compute_velocity,
                      Parameters const& parameters,
                      not_null<std::vector<SystemState>*> const solution) const;

  VanishingCoefficients vanishing_coefficients_;
  // Null if, and only if, |vanishing_coefficients_ == kNone|.
  // If |vanishing_coefficients_ == kFirstBVanishes|, this contains the first
  // and last a coefficients.
  // If |vanishing_coefficients_ == FirstAVanishes|, this contains the first
  // and last b coefficients.
  // These are used to desynchronize and synchronize first-same-as-last
  // integrators.
  std::unique_ptr<FirstSameAsLast> first_same_as_last_;

  // The number of stages, or the number of stages minus one for a
  // first-same-as-last integrator.
  int stages_;

  // The position and momentum nodes.
  // If |vanishing_coefficients_ == kFirstBVanishes|, we do not store the first
  // b coefficient, and the last entry of |a_| is the sum of the first and last
  // a coefficients.
  // If |vanishing_coefficients_ == kLastAVanishes|, we do not store the last a
  // coefficient, and the first entry of |b_| is the sum of the first and last b
  // coefficients.
  std::vector<double> a_;
  std::vector<double> b_;

  // The weights.  Note that the first c coefficient is not stored if
  // |vanishing_coefficients_ == kFirstBVanishes|.
  std::vector<double> c_;
};

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
