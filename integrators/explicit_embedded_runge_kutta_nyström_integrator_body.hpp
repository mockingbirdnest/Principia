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
  // The low-order integration is not performed using compensated summation.
  std::vector<Position> q_low_order(dimension);
  std::vector<Velocity> v_low_order(dimension);
  // TODO(egg): this is a rectangular container, use something more appropriate.
  std::vector<std::vector<Acceleration>> g(stages_);
  for (auto& g_stage : g) {
    g_stage.resize(dimension);
  }
  Time t_stage;
  std::vector<Displacement> q_stage(dimension);

  for (int k = 0; k < dimension; ++k) {
    ∆q[k] *= h;
  }
  for (int i = 0; i < stages_; ++i) {
    t_stage = t.value + c[i] * h;
    for (int k = 0; k < dimension; ++k) {
      for (int j = 0; j < i; ++j) {
        q_stage[k] += g[j][k] * a[i][j]
      }
      q_stage[k] *= h;
      q_stage[k] += c[i] * v[k].value;
      q_stage[k] *= h;
      q_stage[k] += q[k].value;
    }
    compute_acceleration(t_stage, q_stage, &g[i]);
  }
  for (int k = 0; k < dimension; ++k) {
    for (int i = 0; i < stages_; ++i) {
      ∆q[k] += b_hat_[i] * g[i][k];
    }
    ∆q[k] *= h;
    ∆q[k] += v[k].value;
    ∆q[k] *= h;
  }
}

}  // namespace integrators
}  // namespace principia
