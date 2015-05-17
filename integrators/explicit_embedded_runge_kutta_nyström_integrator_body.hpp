#pragma once

#include "integrators/explicit_embedded_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {

using quantities::Difference;
using quantities::Quotient;

namespace integrators {

template<typename Position>
void ExplicitEmbeddedRungeKuttaNyströmIntegrator::Solve(
    RightHandSideComputation<Position> compute_acceleration,
    SystemState<Position, Variation<Position>> const& initial_value,
    Time const& t_final,
    Time const& first_time_step,
    Difference<Position> const& position_tolerance,
    Variation<Position> const& velocity_tolerance,
    double const safety_factor,
    not_null<Solution<Position, Variation<Position>>*> const solution) const {
  int const dimension = initial_value.positions.size();
  CHECK_EQ(initial_value.momenta.size(), dimension);
  CHECK_NE(first_time_step, Time());
  bool const forward = first_time_step > Time();
  if (forward) {
    CHECK_LT(initial_value.time, t_final);
  } else {
    // Integrating backward.
    CHECK_GT(initial_value.time, t_final);
  }
  using Displacement = Difference<Position>;
  using Velocity = Variation<Position>;
  using Acceleration = Variation<Velocity>;
  // The timestep.
  Time h = first_time_step;
  DoublePrecision<Time> t;
  std::vector<Displacement> ∆q(dimension);
  std::vector<Velocity> ∆v(dimension);
  std::vector<DoublePrecision<Position>> q(dimension);
  std::vector<DoublePrecision<Position>> v(dimension);
  std::vector<Position> ∆q_low_order(dimension);
  std::vector<Velocity> ∆v_low_order(dimension);
  // TODO(egg): this is a rectangular container, use something more appropriate.
  std::vector<std::vector<Acceleration>> g(stages_);
  for (auto& g_stage : g) {
    g_stage.resize(dimension);
  }
  Time t_stage;
  std::vector<Position> q_stage(dimension);

  std::fill(∆q.begin(), ∆q.end(), Displacement());

  for (int i = 0; i < stages_; ++i) {
    t_stage = t.value + c[i] * h;
    for (int k = 0; k < dimension; ++k) {
      Acceleration ∑_a_ij_g_j;
      for (int j = 0; j < i; ++j) {
        ∑_a_ij_g_j += g[j][k] * a[i][j]
      }
      q_stage[k] = q[k].value + h * (c[i] * v[k].value + h * ∑_a_ij_g_j);
    }
    compute_acceleration(t_stage, q_stage, &g[i]);
  }

  for (int k = 0; k < dimension; ++k) {
    Acceleration ∑_b_hat_i_g_i;
    Acceleration ∑_b_i_g_i;
    Acceleration ∑_b_prime_hat_i_g_i;
    Acceleration ∑_b_prime_i_g_i;
    // TODO(egg): this would be more readable if it were aligned on the += / =.
    for (int i = 0; i < stages_; ++i) {
      ∑_b_hat_i_g_i += b_hat_[i] * g[i][k];
      ∑_b_i_g_i += b_[i] * g[i][k];
      ∑_b_prime_hat_i_g_i += b_prime_hat_[i] * g[i][k];
      ∑_b_prime_i_g_i += b_prime_[i] * g[i][k];
    }
    ∆q[k] = h * (h * (∑_b_hat_i_g_i) + v[k].value);
    ∆q_low_order[k] = h * (h * (∑_b_i_g_i) + v[k].value);
    ∆v[k] = h * ∑_b_prime_hat_i_g_i;
    ∆v_low_order[k] = h * ∑_b_prime_i_g_i;
  }
}

}  // namespace integrators
}  // namespace principia
