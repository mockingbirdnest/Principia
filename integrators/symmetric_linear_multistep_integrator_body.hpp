#pragma once

#include "integrators/symmetric_linear_multistep_integrator.hpp"

namespace principia {
namespace integrators {

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::
SymmetricLinearMultistepIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind,
    FixedVector<double, half_order_> const & ɑ,
    FixedVector<double, half_order_> const & β_numerators,
    double const β_denominator)
    : FixedStepSizeIntegrator(kind),
      ɑ_(ɑ) {
  for (int i = 0; i < half_order_; ++i) {
    β_[i] = β_numerators[i] / β_denominator;
  }
}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::Solve(
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

  typename ODE::SystemState current_state = *problem.initial_state;

  // Time step.
  Time const& h = step;
  // Current time.
  DoublePrecision<Instant>& t = current_state.time;
  // Current position.
  std::vector<DoublePrecision<Position>>& q = current_state.positions;

  while (h <= Abs((problem.t_final - t.value) - t.error)) {
  }
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8A() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8A,
      {1.0, -2.0, 2.0, -2.0, 2.0},
      {0.0, 22081.0, -29418.0, 75183.0, -75212.0},
      15120.0);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8B() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8B,
      {1.0, 0.0, 0.0, -1 / 2.0, -1.0},
      {0.0, 192481.0, 6582.0, 816783.0, -156812.0},
      120960.0);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const&
QuinlanTremaine1990Order8() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_8,
      {1.0, -2.0, 2.0, -1.0, 0.0},
      {0.0, 17671.0, -23622.0, 61449.0, 50516.0},
      12096.0);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 10> const&
QuinlanTremaine1990Order10() {
  static SymmetricLinearMultistepIntegrator<Position, 10> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_10,
      {1.0, -1.0, 1.0, -1.0, 1.0, 2.0},
      {0.0, 399187.0, -485156.0, 2391436.0, -2816732.0, 4651330.0},
      241920.0);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 12> const&
QuinlanTremaine1990Order12() {
  static SymmetricLinearMultistepIntegrator<Position, 12> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_12,
      {1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0},
      {0.0,
       90987349.0,
       -229596838.0,
       812627169.0,
       -1628539944.0,
       2714971338.0,
       -3041896548.0},
      53222400.0);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 14> const&
QuinlanTremaine1990Order14() {
  static SymmetricLinearMultistepIntegrator<Position, 14> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_14,
      {1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0,
       433489274083.0,
       -1364031998256.0,
       5583113380398.0,
       -14154444148720.0,
       28630585332045.0,
       -42056933842656.0,
       48471792742212.0},
      237758976000.0);
  return integrator;
}

}  // namespace integrators
}  // namespace principia
