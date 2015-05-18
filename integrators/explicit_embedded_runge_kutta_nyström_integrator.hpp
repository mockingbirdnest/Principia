#pragma once

#include <functional>
#include <vector>

#include "base/not_null.hpp"
#include "integrators/motion_integrator.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using quantities::Time;
using quantities::Variation;

namespace integrators {

class ExplicitEmbeddedRungeKuttaNyströmIntegrator {
 public:
  ExplicitEmbeddedRungeKuttaNyströmIntegrator(
      std::vector<double> const& c,
      std::vector<std::vector<double>> const& a,
      std::vector<double> const& b_hat,
      std::vector<double> const& b_prime_hat,
      std::vector<double> const& b,
      std::vector<double> const& b_prime);

  virtual ~ExplicitEmbeddedRungeKuttaNyströmIntegrator() = default;

  ExplicitEmbeddedRungeKuttaNyströmIntegrator() = delete;
  ExplicitEmbeddedRungeKuttaNyströmIntegrator(
      ExplicitEmbeddedRungeKuttaNyströmIntegrator const&) = delete;
  ExplicitEmbeddedRungeKuttaNyströmIntegrator(
      ExplicitEmbeddedRungeKuttaNyströmIntegrator&&) = delete;
  ExplicitEmbeddedRungeKuttaNyströmIntegrator& operator=(
      ExplicitEmbeddedRungeKuttaNyströmIntegrator const&) = delete;
  ExplicitEmbeddedRungeKuttaNyströmIntegrator& operator=(
      ExplicitEmbeddedRungeKuttaNyströmIntegrator&&) = delete;

  // TODO(egg): Copied from MotionIntegrator, unify.
  template<typename Position, typename Momentum>
  struct SystemState {
    std::vector<DoublePrecision<Position>> positions;
    std::vector<DoublePrecision<Momentum>> momenta;
    DoublePrecision<Time> time;
  };

  template<typename Position, typename Momentum>
  using Solution = std::vector<SystemState<Position, Momentum>>;

  template<typename Position>
  using RightHandSideComputation =
      std::function<
          void(Time const& t,
               std::vector<Position> const&,
               not_null<std::vector<Variation<Variation<Position>>>*> const)>;

  template<typename Position>
  using StepSizeController =
      std::function<
          double(std::vector<Difference<Position>> const& q_error_estimate,
                 std::vector<Variation<Position>> const& v_error_estimate)>;

  // TODO(egg): maybe wrap that in some sort of Parameters struct when it's
  // unified with the other integrators.
  template<typename Position>
  void Solve(
      RightHandSideComputation<Position> compute_acceleration,
      SystemState<Position, Variation<Position>> const& initial_value,
      Time const& t_final,
      Time const& first_time_step,
      StepSizeController<Position> step_size_controller,
      double const safety_factor,
      not_null<Solution<Position, Variation<Position>>*> const solution) const;

 protected:
  int stages_;
  int lower_order_;
  // The Runge-Kutta matrix.
  // TODO(egg): This is really a strictly lower-triangular matrix, so we should
  // store it in a smarter way eventually.
  std::vector<std::vector<double>> a_;
  // The nodes.
  std::vector<double> c_;
  // The weights b̂ for the high-order method for the positions.
  std::vector<double> b_hat_;
  // The weights b̂′ for the high-order method for the velocities.
  std::vector<double> b_prime_hat_;
  // The weights b for the low-order method for the positions.
  std::vector<double> b_;
  // The weights b′ for the low-order method for the velocities.
  std::vector<double> b_prime_;
};

// Coefficients form Dormand, El-Mikkawy and Prince (1986),
// Families of Runge-Kutta-Nyström formulae, table 3 (the RK4(3)4FM).
// high order: 4;
// low order:  3;
// stages:     4;
// first-same-as-last, minimizes the 4th order truncation error.
ExplicitEmbeddedRungeKuttaNyströmIntegrator const&
DormandElMikkawyPrince1986RKN434FM();

}  // namespace integrators
}  // namespace principia

#include "integrators/explicit_embedded_runge_kutta_nyström_integrator_body.hpp"
