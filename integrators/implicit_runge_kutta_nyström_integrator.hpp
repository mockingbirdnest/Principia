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
      std::vector<std::vector<double>> const& a,
      std::vector<double> const& b,
      std::vector<double> const& c);

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
  void Solve(RightHandSideComputation<Position> compute_acceleration,
      SystemState<Position, Variation<Position>> const& initial_value,
      Time const& t_final,
      Time const& first_time_step,
      Difference<Position> const& position_tolerance,
      Variation<Position> const& velocity_tolerance,
      not_null<Solution<Position, Variation<Position>>*> const solution) const;

 protected:
  int stages_;
  // The Runge-Kutta matrix.
  std::vector<std::vector<double>> a_;
  // The nodes.
  std::vector<double> c_;
  // The weights for the high-order method for the positions.
  std::vector<double> b_hat_;
  // The weights for the high-order method for the velocities.
  std::vector<double> b_prime_hat_;
  // The weights for the low-order method for the positions.
  std::vector<double> b_;
  // The weights for the low-order method for the velocities.
  std::vector<double> b_prime_;

}

}  // namespace integrators
}  // namespace principia

#include "integrators/implicit_runge_kutta_nyström_integrator_body.hpp"
