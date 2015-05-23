#pragma once

#include <functional>
#include <vector>

#include "base/not_null.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using quantities::Time;
using quantities::Variation;

namespace integrators {

// This class solves ordinary differential equations of the form q″ = f(q, t)
// using an embedded Runge-Kutta-Nyström method.  We follow the standard
// conventions for the coefficients, i.e.,
//   c for the nodes;
//   a for the Runge-Kutta matrix;
//   b̂ for the position weights of the high-order method;
//   b̂′ for the velocity weights of the high-order method;
//   b for the position weights of the low-order method;
//   b′ for the velocity weights of the low-order method.
// See Dormand, El-Mikkawy and Prince (1986),
// Families of Runge-Kutta-Nyström formulae, for an example.
// In the implementation, we follow Dormand, El-Mikkawy and Prince in calling
// the results of the right-hand-side evaluations gᵢ.  The notations kᵢ or fᵢ
// also appear in the litterature.
// Since we are interested in physical applications, we call the solution q and
// its derivative v, rather than the more common y and y′ found in the
// literature on Runge-Kutta-Nyström methods.

class EmbeddedExplicitRungeKuttaNyströmIntegrator {
 public:
  EmbeddedExplicitRungeKuttaNyströmIntegrator(
      std::vector<double> const& c,
      std::vector<std::vector<double>> const& a,
      std::vector<double> const& b_hat,
      std::vector<double> const& b_prime_hat,
      std::vector<double> const& b,
      std::vector<double> const& b_prime,
      int const lower_order);

  ~EmbeddedExplicitRungeKuttaNyströmIntegrator() = default;

  EmbeddedExplicitRungeKuttaNyströmIntegrator() = delete;
  EmbeddedExplicitRungeKuttaNyströmIntegrator(
      EmbeddedExplicitRungeKuttaNyströmIntegrator const&) = delete;
  EmbeddedExplicitRungeKuttaNyströmIntegrator(
      EmbeddedExplicitRungeKuttaNyströmIntegrator&&) = delete;
  EmbeddedExplicitRungeKuttaNyströmIntegrator& operator=(
      EmbeddedExplicitRungeKuttaNyströmIntegrator const&) = delete;
  EmbeddedExplicitRungeKuttaNyströmIntegrator& operator=(
      EmbeddedExplicitRungeKuttaNyströmIntegrator&&) = delete;

  template<typename Position>
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;

  template<typename Position>
  void Solve(
      IntegrationProblem<ODE<Position>> const& problem,
      AdaptiveStepSize<ODE<Position>> const& adaptive_step_size) const;

 protected:
  int stages_;
  int lower_order_;
  std::vector<double> c_;
  // TODO(egg): This is really a strictly lower-triangular matrix, so we should
  // store it in a smarter way eventually.
  std::vector<std::vector<double>> a_;
  std::vector<double> b_hat_;
  std::vector<double> b_prime_hat_;
  std::vector<double> b_;
  std::vector<double> b_prime_;
};

// Coefficients form Dormand, El-Mikkawy and Prince (1986),
// Families of Runge-Kutta-Nyström formulae, table 3 (the RK4(3)4FM).
// high order: 4;
// low order:  3;
// stages:     4;
// first-same-as-last, minimizes the 4th order truncation error.
EmbeddedExplicitRungeKuttaNyströmIntegrator const&
DormandElMikkawyPrince1986RKN434FM();

}  // namespace integrators
}  // namespace principia

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator_body.hpp"
