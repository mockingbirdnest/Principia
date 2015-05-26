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
  // Accelerations at each stage.
  // TODO(egg): this is a rectangular container, use something more appropriate.
  std::vector<std::vector<Acceleration>> g(stages);
  for (auto& g_stage : g) {
    g_stage.resize(dimension);
  }

  bool at_end = false;
}

}  // namespace integrators
}  // namespace principia
