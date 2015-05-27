#pragma once

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {

template<typename Position, int order, int evaluations,
         CompositionMethod composition>
SymplecticRungeKuttaNyströmIntegrator<Position, order, evaluations,
                                      composition>::
SymplecticRungeKuttaNyströmIntegrator(FixedVector<double, stages_> const& a,
                                      FixedVector<double, stages_> const& b,
                                      FixedVector<double, stages_> const& c)
    : a_(a), b_(b), c_(c) {}

template<typename Position, int order, int evaluations,
         CompositionMethod composition>
void SymplecticRungeKuttaNyströmIntegrator<Position, order, evaluations,
                                           composition>::Solve(
    IntegrationProblem<ODE> const& problem,
    Time const& step) const {
  using Displacement = typename ODE::Displacement;
  using Velocity = typename ODE::Velocity;
  using Acceleration = typename ODE::Acceleration;

  // Argument checks.
  CHECK_NOTNULL(problem.initial_state);
  int const dimension = problem.initial_state->positions.size();
  CHECK_EQ(dimension, problem.initial_state->velocities.size());
  CHECK_NE(Time(), step);
  Sign const integration_direction = Sign(step);
  if (integration_direction.Positive()) {
    // Integrating forward.
    CHECK_LT(problem.initial_state->time.value, problem.t_final);
  } else {
    // Integrating backward.
    CHECK_GT(problem.initial_state->time.value, problem.t_final);
  }

  typename ODE::SystemState current_state = *problem.initial_state;

  // Time step.
  Time h = step;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<Instant>& t = current_state.time;

  // Position increment.
  std::vector<Displacement> ∆q(dimension);
  // Velocity increment.
  std::vector<Velocity> ∆v(dimension);
  // Current position.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Position>>& q = current_state.positions;
  // Current velocity.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Velocity>>& v = current_state.velocities;

  // Current Runge-Kutta-Nyström stage.
  std::vector<Position> q_stage(dimension);
  // Accelerations at the current stage.
  std::vector<Acceleration> g(dimension);

  bool at_end = false;

  for (;;) {
    // TODO(egg): do we want to keep cover the whole integration interval, or
    // should the last point be in the interval?  This implements the latter.
    // There's also the possibility of making the caller determine the
    // termination condition.
    Time const time_to_end = (problem.t_final - t.value) - t.error;
    at_end = integration_direction * h > integration_direction * time_to_end;
    if (at_end) {
      break;
    }

    std::fill(∆q.begin(), ∆q.end(), Displacement{});
    std::fill(∆v.begin(), ∆v.end(), Velocity{});
    for (int i = 0; i < stages_; ++i) {
      for (int k = 0; k < dimension; ++k) {
        q_stage[k] = q[k].value + ∆q[k];
      }
      problem.equation.compute_acceleration(q_stage, t.value + c_[i] * h, &g);
      for (int k = 0; k < dimension; ++k) {
        // TODO(egg): reformulate to reduce roundoff error.
        ∆v[k] += h * b_[i] * g[k];
        ∆q[k] += h * a_[i] * (v[k].value + ∆v[k]);
      }
    }

    // Increment the solution.
    t.Increment(h);
    for (int k = 0; k < dimension; ++k) {
      q[k].Increment(∆q[k]);
      v[k].Increment(∆v[k]);
    }
    problem.append_state(current_state);
  }

}

}  // namespace integrators
}  // namespace principia
