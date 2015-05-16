#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

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

// We follow the convention of McLachlan and Atela, calling the position nodes
// aᵢ and the momentum nodes bᵢ.  The following notations appear in the
// litterature:
//   (a, b) in McLachlan and Atela, as well as Candy and Rozmus;
//   (d, c) in Ruth, Yoshida, as well as Forest and Ruth;
//   (B, b) in Sofroniou and Spaletta.
// Moreover, we follow the convention of Sofroniou and Spaletta in calling the
// weights used for the time argument of the force computation cᵢ, with
// c₁ = 0, cᵢ = cᵢ₋₁ + aᵢ₋₁ for i > 1.

class SPRKIntegrator : public SRKNIntegrator {
 public:
  // NOTE(egg): this should be a using when we have VS 2015.
  SPRKIntegrator(std::vector<double> const& a, std::vector<double> const& b);

  virtual ~SPRKIntegrator() = default;

  SPRKIntegrator() = delete;
  SPRKIntegrator(SPRKIntegrator const&) = delete;
  SPRKIntegrator(SPRKIntegrator&&) = delete;
  SPRKIntegrator& operator=(SPRKIntegrator const&) = delete;
  SPRKIntegrator& operator=(SPRKIntegrator&&) = delete;

  template<typename Position, typename Momentum>
  using SPRKRightHandSideComputation =
      std::function<void(Time const&,
                         std::vector<Position> const&,
                         std::vector<Variation<Momentum>>* const)>;

  template<typename Position, typename Momentum>
  using SPRKAutonomousRightHandSideComputation =
      std::function<void(std::vector<Momentum> const&,
                         std::vector<Variation<Position>>* const)>;

  // The functors |compute_velocity| and |compute_force| compute
  // ∂T/∂pᵢ(p) and ∂V/∂qᵢ(q,t) respectively.
  template<typename Position, typename Momentum>
  void SolveIncrement(
      SPRKRightHandSideComputation<Position, Momentum> compute_force,
      SPRKAutonomousRightHandSideComputation<Position, Momentum>
          compute_velocity,
      Parameters<Position, Momentum> const& parameters,
      not_null<Solution<Position, Momentum>*> const solution) const;

 private:
  template<VanishingCoefficients vanishing_coefficients,
           typename Position, typename Momentum>
  void SolveIncrementOptimized(
      SPRKRightHandSideComputation<Position, Momentum> compute_force,
      SPRKAutonomousRightHandSideComputation<Position, Momentum>
          compute_velocity,
      Parameters<Position, Momentum> const& parameters,
      not_null<Solution<Position, Momentum>*> const solution) const;
};

// First-same-as-last (FSAL) methods with k stages require k - 1 force
// and velocity evaluations for sparse output.  For dense output, methods
// with synchronous momenta require k - 1 force evaluations and k velocity
// evaluations, methods with synchronous positions require k force evaluations
// and k - 1 velocity evaluations.

// Second order, 2 stages, FSAL (synchronous momenta).
SPRKIntegrator const& Leapfrog();
// Second order, 2 stages, FSAL (synchronous positions).
SPRKIntegrator const& PseudoLeapfrog();
// Second order, 2 stages.  This method minimizes the error constant.
// Coefficients from McLachlan and Atela (1992),
// The accuracy of symplectic integrators, table 2.
// http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
SPRKIntegrator const& McLachlanAtela1992Order2Optimal();
// Second order, 3 stages, FSAL (synchronous momenta).
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
SPRKIntegrator const& McLachlan1995S2();
// Third order, 3 stages.
// Coefficients from Ruth (1983), A canonical integration technique,
// https://accelconf.web.cern.ch/accelconf/p83/PDF/PAC1983_2669.PDF.
SPRKIntegrator const& Ruth1983();
// Third order, 3 stages.  This method minimizes the error constant.
// Coefficients from McLachlan and Atela (1992),
// The accuracy of symplectic integrators, equation 3.3.
// http://eaton.math.rpi.edu/CSUMS/Papers/Symplectic/McLachlan_Atela_92.pdf.
// NOTE(egg): the coefficients given in table 2 for this integrator are
// incorrect (b1 and b2 are swapped).
SPRKIntegrator const& McLachlanAtela1992Order3Optimal();
// Fourth order, 4 stages, FSAL (synchronous momenta).
// Coefficients from Forest and Ruth (1990),
// Fourth-order symplectic integration, equation 4.8.
// http://zwe.web.cern.ch/zwe/CAS/biblio/ruth-forest.pdf.
// This scheme was independently discovered by Candy and Rozmus (1991),
// A Symplectic Integration Algorithm for Separable Hamiltonian Functions
// (submitted earlier and published later than the Forest and Ruth paper).
SPRKIntegrator const& CandyRozmus1991ForestRuth1990SynchronousMomenta();
// Fourth order, 4 stages, FSAL (synchronous positions).
// Coefficients from Forest and Ruth (1990),
// Fourth-order symplectic integration, equation 4.9.
// http://zwe.web.cern.ch/zwe/CAS/biblio/ruth-forest.pdf.
// This scheme was independently discovered by Candy and Rozmus (1991),
// A Symplectic Integration Algorithm for Separable Hamiltonian Functions
// (submitted earlier and published later than the Forest and Ruth paper).
SPRKIntegrator const& CandyRozmus1991ForestRuth1990SynchronousPositions();
// Fourth order, 6 stages, FSAL (synchronous momenta).
// Coefficients from Suzuki (1990), Fractal decompositions of exponential
// operators with applications to many-body theories and Monte Carlo
// simulations,
SPRKIntegrator const& Suzuki1990();
// Fourth order, 6 stages, FSAL (synchronous momenta).
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
SPRKIntegrator const& McLachlan1995SS5();
// Fourth order, 5 stages, FSAL (synchronous momenta).  Ibidem.
SPRKIntegrator const& McLachlan1995S4();
// Fourth order, 6 stages, FSAL (synchronous momenta).  Ibidem.
SPRKIntegrator const& McLachlan1995S5();
// Fourth order, 7 stages, FSAL (synchronous momenta).
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
SPRKIntegrator const& BlanesMoan2002S6();
// Sixth order, 8 stages, FSAL (synchronous positions).
// Coefficients from Yoshida (1990),
// Construction of higher order symplectic integrators
// http://sixtrack.web.cern.ch/SixTrack/doc/yoshida00.pdf.
// NOTE(egg): The coefficients were derived from equations 5.4 through 5.17
// rather than computed from the wᵢ given in tables 1 and 2.  The results were
// then cross-checked against those obtained from the tables.
SPRKIntegrator const& Yoshida1990Order6A();
// Sixth order, 8 stages, FSAL (synchronous positions).  Ibidem.
SPRKIntegrator const& Yoshida1990Order6B();
// Sixth order, 8 stages, FSAL (synchronous positions).  Ibidem.
SPRKIntegrator const& Yoshida1990Order6C();
// Sixth order, 10 stages, FSAL (synchronous momenta).
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
SPRKIntegrator const& McLachlan1995SS9();
// Sixth order, 11 stages, FSAL (synchronous momenta).
// Coefficients from Blanes and Moan (2002),
// Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods,
// http://personales.upv.es/serblaza/2002JCAM.pdf.
SPRKIntegrator const& BlanesMoan2002S10();
// Eighth order, 16 stages, FSAL (synchronous positions).
// Coefficients from Yoshida (1990),
// Construction of higher order symplectic integrators
// http://sixtrack.web.cern.ch/SixTrack/doc/yoshida00.pdf.
SPRKIntegrator const& Yoshida1990Order8A();
// Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
SPRKIntegrator const& Yoshida1990Order8B();
// Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
SPRKIntegrator const& Yoshida1990Order8C();
// Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
SPRKIntegrator const& Yoshida1990Order8D();
// Eighth order, 16 stages, FSAL (synchronous positions).  Ibidem.
SPRKIntegrator const& Yoshida1990Order8E();
// Eighth order, 16 stages, FSAL (synchronous momenta).
// Coefficients from McLachlan (1995),
// On the numerical integration of ordinary differential equations by symmetric
// composition methods, http://www.massey.ac.nz/~rmclachl/sisc95.pdf.
SPRKIntegrator const& McLachlan1995SS15();
// Eighth order, 18 stages, FSAL (synchronous momenta).  Ibidem.
SPRKIntegrator const& McLachlan1995SS17();

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_partitioned_runge_kutta_integrator_body.hpp"
