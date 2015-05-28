#pragma once

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"

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

template<typename Position>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Position, 4, 3, 4, true> const&
DormandElMikkawyPrince1986RKN434FM() {
  static EmbeddedExplicitRungeKuttaNyströmIntegrator<
             Position, 4, 3, 4, true> const integrator(
      // c
      { 0.0         ,   1.0 /   4.0,   7.0 /  10.0,  1.0},
      // a
      {
        1.0 /   32.0,
        7.0 / 1000.0, 119.0 / 500.0,
        1.0 /   14.0,   8.0 /  27.0,  25.0 / 189.0},
      // b̂
      { 1.0 /   14.0,   8.0 /  27.0,  25.0 / 189.0,  0.0},
      // b̂′
      { 1.0 /   14.0,  32.0 /  81.0, 250.0 / 567.0,  5.0 / 54.0},
      // b
      {-7.0 /  150.0,  67.0 / 150.0,   3.0 /  20.0, -1.0 / 20.0},
      // b′
      {13.0 /   21.0, -20.0 /  27.0, 275.0 / 189.0, -1.0 /  3.0});
  return integrator;
}

template<typename Position, int higher_order, int lower_order, int stages,
         bool first_same_as_last>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Position, higher_order, lower_order,
                                            stages, first_same_as_last>::
EmbeddedExplicitRungeKuttaNyströmIntegrator(
    FixedVector<double, stages> const& c,
    FixedStrictlyLowerTriangularMatrix<double, stages> const& a,
    FixedVector<double, stages> const& b_hat,
    FixedVector<double, stages> const& b_prime_hat,
    FixedVector<double, stages> const& b,
    FixedVector<double, stages> const& b_prime)
    : c_(c),
      a_(a),
      b_hat_(b_hat),
      b_prime_hat_(b_prime_hat),
      b_(b),
      b_prime_(b_prime) {
  // the first node is always 0 in an explicit method.
  CHECK_EQ(0.0, c_[0]);
  if (first_same_as_last) {
    // Check that the conditions for the FSAL property are satisfied, see for
    // instance Dormand, El-Mikkawy and Prince (1986),
    // Families of Runge-Kutta-Nyström formulae, equation 3.1.
    CHECK_EQ(1.0, c_[stages - 1]);
    CHECK_EQ(0.0, b_hat_[stages - 1]);
    for (int j = 0; j < stages - 1; ++j) {
      CHECK_EQ(b_hat_[j], a_[stages - 1][j]);
    }
  }
}

template<typename Position, int higher_order, int lower_order, int stages,
         bool first_same_as_last>
void EmbeddedExplicitRungeKuttaNyströmIntegrator<Position,
                                                 higher_order,
                                                 lower_order,
                                                 stages,
                                                 first_same_as_last>::Solve(
    IntegrationProblem<ODE> const& problem,
    AdaptiveStepSize<ODE> const& adaptive_step_size) const {
  using Displacement = typename ODE::Displacement;
  using Velocity = typename ODE::Velocity;
  using Acceleration = typename ODE::Acceleration;

  // Argument checks.
  CHECK_NOTNULL(problem.initial_state);
  int const dimension = problem.initial_state->positions.size();
  CHECK_EQ(dimension, problem.initial_state->velocities.size());
  CHECK_NE(Time(), adaptive_step_size.first_time_step);
  Sign const integration_direction =
      Sign(adaptive_step_size.first_time_step);
  if (integration_direction.Positive()) {
    // Integrating forward.
    CHECK_LT(problem.initial_state->time.value, problem.t_final);
  } else {
    // Integrating backward.
    CHECK_GT(problem.initial_state->time.value, problem.t_final);
  }
  CHECK_GT(adaptive_step_size.safety_factor, 0);
  CHECK_LT(adaptive_step_size.safety_factor, 1);

  typename ODE::SystemState current_state = *problem.initial_state;

  // Time step.
  Time h = adaptive_step_size.first_time_step;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<Instant>& t = current_state.time;

  // Position increment (high-order).
  std::vector<Displacement> ∆q_hat(dimension);
  // Velocity increment (high-order).
  std::vector<Velocity> ∆v_hat(dimension);
  // Current position.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Position>>& q_hat = current_state.positions;
  // Current velocity.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Velocity>>& v_hat = current_state.velocities;

  // Difference between the low- and high-order approximations.
  typename ODE::SystemStateError error_estimate;
  error_estimate.position_error.resize(dimension);
  error_estimate.velocity_error.resize(dimension);

  // Current Runge-Kutta-Nyström stage.
  std::vector<Position> q_stage(dimension);
  // Accelerations at each stage.
  // TODO(egg): this is a rectangular container, use something more appropriate.
  std::vector<std::vector<Acceleration>> g(stages);
  for (auto& g_stage : g) {
    g_stage.resize(dimension);
  }

  bool at_end = false;
  double tolerance_to_error_ratio;

  // The first stage of the Runge-Kutta-Nyström iteration.  In the FSAL case,
  // |first_stage == 1| after the first step, since the first RHS evaluation has
  // already occured in the previous step.  In the non-FSAL case and in the
  // first step of the FSAL case, |first_stage == 0|.
  int first_stage = 0;

  // No step size control on the first step.
  goto runge_kutta_nyström_step;

  while (!at_end) {
    // Compute the next step with decreasing step sizes until the error is
    // tolerable.
    do {
      // Adapt step size.
      // TODO(egg): find out whether there's a smarter way to compute that root,
      // especially if we make the order compile-time.
      h *= adaptive_step_size.safety_factor *
               std::pow(tolerance_to_error_ratio, 1.0 / (lower_order + 1));

    runge_kutta_nyström_step:
      // Termination condition.
      Time const time_to_end = (problem.t_final - t.value) - t.error;
      at_end = integration_direction * h >= integration_direction * time_to_end;
      if (at_end) {
        // The chosen step size will overshoot.  Clip it to just reach the end,
        // and terminate if the step is accepted.
        h = time_to_end;
      }

      // Runge-Kutta-Nyström iteration; fills |g|.
      for (int i = first_stage; i < stages; ++i) {
        Instant const t_stage = t.value + c_[i] * h;
        for (int k = 0; k < dimension; ++k) {
          Acceleration ∑j_a_ij_g_jk{};
          for (int j = 0; j < i; ++j) {
            ∑j_a_ij_g_jk += a_[i][j] * g[j][k];
          }
          q_stage[k] = q_hat[k].value +
                           h * (c_[i] * v_hat[k].value + h * ∑j_a_ij_g_jk);
        }
        problem.equation.compute_acceleration(t_stage, q_stage, &g[i]);
      }

      // Increment computation and step size control.
      for (int k = 0; k < dimension; ++k) {
        Acceleration ∑i_b_hat_i_g_ik{};
        Acceleration ∑i_b_i_g_ik{};
        Acceleration ∑i_b_prime_hat_i_g_ik{};
        Acceleration ∑i_b_prime_i_g_ik{};
        // Please keep the eight assigments below aligned, they become illegible
        // otherwise.
        for (int i = 0; i < stages; ++i) {
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

        error_estimate.position_error[k] = ∆q_k - ∆q_hat[k];
        error_estimate.velocity_error[k] = ∆v_k - ∆v_hat[k];
      }
      tolerance_to_error_ratio =
          adaptive_step_size.tolerance_to_error_ratio(h, error_estimate);
    } while (tolerance_to_error_ratio < 1.0);

    if (first_same_as_last) {
      using std::swap;
      swap(g.front(), g.back());
      first_stage = 1;
    }

    // Increment the solution with the high-order approximation.
    t.Increment(h);
    for (int k = 0; k < dimension; ++k) {
      q_hat[k].Increment(∆q_hat[k]);
      v_hat[k].Increment(∆v_hat[k]);
    }
    problem.append_state(current_state);
  }
}

}  // namespace integrators
}  // namespace principia
