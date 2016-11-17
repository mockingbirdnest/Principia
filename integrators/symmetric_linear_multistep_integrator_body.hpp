#pragma once

#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::
SymmetricLinearMultistepIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind,
    FixedStepSizeIntegrator<ODE> const& startup_integrator,
    FixedVector<double, half_order_> const & ɑ,
    FixedVector<double, half_order_> const& β_numerators,
    double const β_denominator)
    : FixedStepSizeIntegrator(kind),
      startup_integrator_(startup_integrator),
      ɑ_(ɑ),
      β_numerators_(β_numerators),
      β_denominator_(β_denominator) {}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::Solve(
    Instant const& t_final,
    not_null<IntegrationInstance*> const instance) const {
  using Displacement = typename ODE::Displacement;
  using Velocity = typename ODE::Velocity;
  using Acceleration = typename ODE::Acceleration;

  Instance* const down_cast_instance =
      dynamic_cast_not_null<Instance*>(instance);
  auto const& equation = down_cast_instance->equation;
  auto const& append_state = down_cast_instance->append_state;
  Time const& step = down_cast_instance->step;

  CHECK(!current_state.empty());
  if (current_states.size() < order_ - 1) {
    auto const startup_append_state =
        [&current_states](typename ODE::SystemState const& state) {
      current_states.push_back(state);
      if (current_states.size() == order_ - 1) {
        current_states.pop_front();
      }
    };
    typename ODE::SystemState const& startup_initial_state =
        current_states.back();
    auto const startup_instance =
        startup_integrator_.NewInstance({equation, &startup_initial_state},
                                        startup_append_state,
                                        step);
    startup_integrator_.Solve(
        std::min(current_states.back().time.value +
                     (order_ - 1 - current_states.size()) * step,
                 t_final),
        startup_instance.get());
  }
  CHECK_EQ(current_states.size(), order_ - 1);

  // Argument checks.
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

  while (h <= Abs((t_final - t.value) - t.error)) {
  }
}

template<typename Position, int order_>
not_null<std::unique_ptr<IntegrationInstance>>
SymmetricLinearMultistepIntegrator<Position, order_>::NewInstance(
    IntegrationProblem<ODE> const& problem,
    IntegrationInstance::AppendState<ODE> append_state,
    Time const& step) const {}

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::Instance::Instance(
    IntegrationProblem<ODE> problem,
    AppendState<ODE> append_state,
    Time step) {}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::StartupSolve(
    Instant const& t_final,
    Instance& instance) const {
  CHECK(!instance.current_states.empty());
  CHECK_LT(instance.current_states.size(), order_ - 1);

  auto const startup_append_state =
      [&instance](typename ODE::SystemState const& state) {
        instance.current_states.push_back(state);
        if (instance.current_states.size() == order_ - 1) {
          instance.current_states.pop_front();
        }
      };
  typename ODE::SystemState const& startup_initial_state =
      instance.current_states.back();
  auto const startup_instance = startup_integrator_.NewInstance(
      {instance.equation, &startup_initial_state},
      startup_append_state,
      instance.step);

  startup_integrator_.Solve(
      std::min(
          instance.current_states.back().time.value +
              (order_ - 1 - instance.current_states.size()) * instance.step,
          t_final),
      startup_instance.get());

  CHECK_EQ(current_states.size(), order_ - 1);
}

// TODO(phl): Pick more appropriate startup integrators below.
template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8A() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8A,
      McLachlanAtela1992Order5Optimal<Position>(),
      {1.0, -2.0, 2.0, -2.0, 2.0},
      {0.0, 22081.0, -29418.0, 75183.0, -75212.0},
      15120.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8B() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8B,
      McLachlanAtela1992Order5Optimal<Position>(),
      {1.0, 0.0, 0.0, -1 / 2.0, -1.0},
      {0.0, 192481.0, 6582.0, 816783.0, -156812.0},
      120960.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const&
QuinlanTremaine1990Order8() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_8,
      McLachlanAtela1992Order5Optimal<Position>(),
      {1.0, -2.0, 2.0, -1.0, 0.0},
      {0.0, 17671.0, -23622.0, 61449.0, 50516.0},
      12096.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 10> const&
QuinlanTremaine1990Order10() {
  static SymmetricLinearMultistepIntegrator<Position, 10> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_10,
      McLachlanAtela1992Order5Optimal<Position>(),
      {1.0, -1.0, 1.0, -1.0, 1.0, 2.0},
      {0.0, 399187.0, -485156.0, 2391436.0, -2816732.0, 4651330.0},
      241920.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 12> const&
QuinlanTremaine1990Order12() {
  static SymmetricLinearMultistepIntegrator<Position, 12> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_12,
      McLachlanAtela1992Order5Optimal<Position>(),
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

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 14> const&
QuinlanTremaine1990Order14() {
  static SymmetricLinearMultistepIntegrator<Position, 14> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_14,
      McLachlanAtela1992Order5Optimal<Position>(),
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
