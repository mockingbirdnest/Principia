#pragma once

#include "integrators/explicit_embedded_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {

using geometry::Sign;
using quantities::Difference;
using quantities::Quotient;

namespace integrators {

inline ExplicitEmbeddedRungeKuttaNyströmIntegrator const&
DormandElMikkawyPrince1986RKN434FM() {
  static ExplicitEmbeddedRungeKuttaNyströmIntegrator const integrator(
      // c
      {0.0, 1.0 / 4.0, 7.0 / 10.0, 1.0},
      // a
      {{},
       {1.0 /   32.0},
       {7.0 / 1000.0, 119.0 / 500.0},
       {1.0 /   14.0,   8.0 /  27.0, 25.0 / 189.0}},
      // b̂
      { 1.0 /  14.0,   8.0 /  27.0,  25.0 / 189.0,  0.0},
      // b̂′
      { 1.0 /  14.0,  32.0 /  81.0, 250.0 / 567.0,  5.0 / 54.0},
      // b
      {-7.0 / 150.0,  67.0 / 150.0,   3.0 /  20.0, -1.0 / 20.0},
      // b′
      {13.0 /  21.0, -20.0 /  27.0, 275.0 / 189.0, -1.0 /  3.0},
      3 /*lower_order*/);
  return integrator;
}

inline ExplicitEmbeddedRungeKuttaNyströmIntegrator::
           ExplicitEmbeddedRungeKuttaNyströmIntegrator(
    std::vector<double> const& c,
    std::vector<std::vector<double>> const& a,
    std::vector<double> const& b_hat,
    std::vector<double> const& b_prime_hat,
    std::vector<double> const& b,
    std::vector<double> const& b_prime,
    int const lower_order)
    : stages_(c.size()),
      lower_order_(lower_order),
      c_(c),
      a_(a),
      b_hat_(b_hat),
      b_prime_hat_(b_prime_hat),
      b_(b),
      b_prime_(b_prime) {
  CHECK_EQ(a_.size(), stages_);
  for(int i = 0; i < stages_; ++i) {
    CHECK_EQ(a_[i].size(), i);
  }
  CHECK_EQ(b_hat_.size(), stages_);
  CHECK_EQ(b_prime_hat_.size(), stages_);
  CHECK_EQ(b_.size(), stages_);
  CHECK_EQ(b_prime_.size(), stages_);
}

template<typename Position>
void ExplicitEmbeddedRungeKuttaNyströmIntegrator::Solve(
    RightHandSideComputation<Position> compute_acceleration,
    SystemState<Position, Variation<Position>> const& initial_value,
    Time const& t_final,
    Time const& first_time_step,
    StepSizeController<Position> step_size_controller,
    double const safety_factor,
    not_null<Solution<Position, Variation<Position>>*> const solution) const {
  using Displacement = Difference<Position>;
  using Velocity = Variation<Position>;
  using Acceleration = Variation<Velocity>;

  // Argument checks.
  int const dimension = initial_value.positions.size();
  CHECK_EQ(initial_value.momenta.size(), dimension);
  CHECK_NE(first_time_step, Time());
  Sign const integration_direction = Sign(first_time_step);
  if (integration_direction.Positive()) {
    // Integrating forward.
    CHECK_LT(initial_value.time.value, t_final);
  } else {
    // Integrating backward.
    CHECK_GT(initial_value.time.value, t_final);
  }
  CHECK_GT(safety_factor, 0);
  CHECK_LT(safety_factor, 1);

  // Time step.
  Time h = first_time_step;
  // Current time.
  DoublePrecision<Time> t = initial_value.time;

  // Position increment (high-order).
  std::vector<Displacement> ∆q_hat(dimension);
  // Velocity increment (high-order).
  std::vector<Velocity> ∆v_hat(dimension);
  // Current position.
  std::vector<DoublePrecision<Position>> q_hat = initial_value.positions;
  // Current velocity.
  std::vector<DoublePrecision<Velocity>> v_hat = initial_value.momenta;

  // Difference between the low- and high-order approximations of the position.
  std::vector<Displacement> q_error_estimate(dimension);
  // Difference between the low- and high-order approximations of the velocity.
  std::vector<Velocity> v_error_estimate(dimension);

  // Current Runge-Kutta-Nyström stage.
  std::vector<Position> q_stage(dimension);
  // Accelerations at each stage.
  // TODO(egg): this is a rectangular container, use something more appropriate.
  std::vector<std::vector<Acceleration>> g(stages_);
  for (auto& g_stage : g) {
    g_stage.resize(dimension);
  }

  bool at_end = false;
  double tolerance_to_error_ratio;

  // No step size control on the first step.
  goto runge_kutta_nyström_step;

  while (!at_end) {
    // Compute the next step with decreasing step sizes until the error is
    // tolerable.
    do {
      // Adapt step size.
      // TODO(egg): find out whether there's a smarter way to compute that root,
      // especially if we make the order compile-time.
      h *= safety_factor * std::pow(tolerance_to_error_ratio,
                                    1.0 / (lower_order_ + 1));

    runge_kutta_nyström_step:
      // Termination condition.
      Time const time_to_end = (t_final - t.value) - t.error;
      at_end = integration_direction * h >= integration_direction * time_to_end;
      if (at_end) {
        // The chosen step size will overshoot.  Clip it to just reach the end,
        // and terminate if the step is accepted.
        h = time_to_end;
      }

      // Runge-Kutta-Nyström iteration; fills |g|.
      for (int i = 0; i < stages_; ++i) {
        Time const t_stage = t.value + c_[i] * h;
        for (int k = 0; k < dimension; ++k) {
          Acceleration ∑j_a_ij_g_jk{};
          for (int j = 0; j < i; ++j) {
            ∑j_a_ij_g_jk += a_[i][j] * g[j][k];
          }
          q_stage[k] = q_hat[k].value +
                           h * (c_[i] * v_hat[k].value + h * ∑j_a_ij_g_jk);
        }
        compute_acceleration(t_stage, q_stage, &g[i]);
      }
      // TODO(egg): handle the FSAL case.

      // Increment computation and step size control.
      for (int k = 0; k < dimension; ++k) {
        Acceleration ∑i_b_hat_i_g_ik{};
        Acceleration ∑i_b_i_g_ik{};
        Acceleration ∑i_b_prime_hat_i_g_ik{};
        Acceleration ∑i_b_prime_i_g_ik{};
        // Please keep the eight assigments below aligned, they become illegible
        // otherwise.
        for (int i = 0; i < stages_; ++i) {
          ∑i_b_hat_i_g_ik       += b_hat_[i] * g[i][k];
          ∑i_b_i_g_ik           += b_[i] * g[i][k];
          ∑i_b_prime_hat_i_g_ik += b_prime_hat_[i] * g[i][k];
          ∑i_b_prime_i_g_ik     += b_prime_[i] * g[i][k];
        }
        // The hat-less ∆q and ∆v are the low-order increments.
        ∆q_hat[k]               = h * (h * (∑i_b_hat_i_g_ik) + v_hat[k].value);
        Displacement const ∆q_k = h * (h * (∑i_b_i_g_ik) + v_hat[k].value);
        ∆v_hat[k]               = h * ∑i_b_prime_hat_i_g_ik;
        Velocity const ∆v_k     = h * ∑i_b_prime_i_g_ik;

        q_error_estimate[k] = ∆q_k - ∆q_hat[k];
        v_error_estimate[k] = ∆v_k - ∆v_hat[k];
      }
      tolerance_to_error_ratio =
          step_size_controller(h, q_error_estimate, v_error_estimate);
    } while (tolerance_to_error_ratio < 1.0);

    // Increment the solution with the high-order approximation.
    t.Increment(h);
    for(int k = 0; k < dimension; ++k) {
      q_hat[k].Increment(∆q_hat[k]);
      v_hat[k].Increment(∆v_hat[k]);
    }
    solution->push_back({q_hat, v_hat, t});
  }
}

}  // namespace integrators
}  // namespace principia
