#pragma once

#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {

namespace {

int const startup_step_divisor = 16;

}  // namespace

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::
SymmetricLinearMultistepIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind,
    FixedStepSizeIntegrator<ODE> const& startup_integrator,
    FixedVector<double, half_order_> const & ɑ,
    FixedVector<double, half_order_> const& β_numerator,
    double const β_denominator)
    : FixedStepSizeIntegrator(kind),
      startup_integrator_(startup_integrator),
      ɑ_(ɑ),
      β_numerator_(β_numerator),
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

  auto& previous_steps = down_cast_instance->previous_steps;

  if (previous_steps.size() < order_ - 1) {
    StartupSolve(t_final, *down_cast_instance);
  }

  // Argument checks.
  int const dimension = previous_steps.back().displacements.size();
  CHECK_LT(Time(), step);

  // Time step.
  Time const& h = step;
  // Current time.
  DoublePrecision<Instant> t = previous_steps.back().time;
  // Order.
  int const k = order_;

  std::vector<DoublePrecision<Position>> Σj_minus_ɑj_qj(dimension);
  std::vector<Acceleration> Σj_βj_numerator_aj(dimension);
  while (h <= (t_final - t.value) - t.error) {
    // We take advantage of the symmetry to iterate on the list of previous
    // steps from both ends.
    auto front_it = previous_steps.begin();
    auto back_it = previous_steps.rbegin();

    // This block corresponds to j = 0.  We must not pair it with j = k.
    {
      std::vector<Displacement> const& qj = front_it->displacements;
      std::vector<Acceleration> const& aj = front_it->accelerations;
      double const ɑj = ɑ_[0];
      double const βj_numerator = β_numerator_[0];
      for (int d = 0; d < dimension; ++d) {
        Σj_minus_ɑj_qj[d] = -ɑj * qj[d];
        Σj_βj_numerator_aj[d] = βj_numerator * aj[d];
      }
      ++front_it;
    }
    // The generic value of j, paired with k - j.
    for (int j = 1; j < k / 2; ++j) {
      std::vector<Displacement> const& qj = front_it->displacements;
      std::vector<Displacement> const& qk_minus_j = back_it->displacements;
      std::vector<Acceleration> const& aj = front_it->accelerations;
      std::vector<Acceleration> const& ak_minus_j = back_it->accelerations;
      double const ɑj = ɑ_[j];
      double const βj_numerator = β_numerator_[j];
      for (int d = 0; d < dimension; ++d) {
        Σj_minus_ɑj_qj[d].Increment(-ɑj * (qj[d] + qk_minus_j[d]));
        Σj_βj_numerator_aj[d] += βj_numerator * (aj[d] + ak_minus_j[d]);
      }
      ++front_it;
      ++back_it;
    }
    // This block corresponds to j = k / 2.  We must not pair it with j = k / 2.
    {
      std::vector<Displacement> const& qj = front_it->displacements;
      std::vector<Acceleration> const& aj = front_it->accelerations;
      double const ɑj = ɑ_[k / 2];
      double const βj_numerator = β_numerator_[k / 2];
      for (int d = 0; d < dimension; ++d) {
        Σj_minus_ɑj_qj[d].Increment(-ɑj * qj[d]);
        Σj_βj_numerator_aj[d] += βj_numerator * aj[d];
      }
    }

    // Create a new step in the instance.
    t.Increment(h);
    previous_steps.pop_front();
    previous_steps.emplace_back();
    Step& current_step = previous_steps.back();
    current_step.time = t;
    current_step.accelerations.resize(dimension);

    // Fill the new step.  We skip the division by ɑk as it is equal to 1.0.
    double const ɑk = ɑ_[0];
    CHECK_EQ(ɑk, 1.0);
    typename ODE::SystemState system_state;
    std::vector<Position> positions;
    positions.reserve(dimension);
    for (int d = 0; d < dimension; ++d) {
      DoublePrecision<Position>& current_position = Σj_minus_ɑj_qj[d];
      current_position.Increment(
          h * h * Σj_βj_numerator_aj[d] / β_denominator_);
      current_step.displacements.push_back(current_position.value - Position());
      positions.push_back(current_position.value);
      //TODO(phl):Improve the computation of the velocities
      auto it = previous_steps.rbegin();
      ++it;
      system_state.velocities.push_back(
          (current_position.value - (Position() + it->displacements[d])) / h);
      system_state.positions.push_back(current_position);
    }
    equation.compute_acceleration(
        t.value, positions, &current_step.accelerations);

    // Inform the caller of the new state.
    system_state.time = t;
    append_state(system_state);
  }
}

template<typename Position, int order_>
not_null<std::unique_ptr<IntegrationInstance>>
SymmetricLinearMultistepIntegrator<Position, order_>::NewInstance(
    IntegrationProblem<ODE> const& problem,
    IntegrationInstance::AppendState<ODE> append_state,
    Time const& step) const {
  return make_not_null_unique<Instance>(problem,
                                        std::move(append_state),
                                        step);
}

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::Instance::Instance(
    IntegrationProblem<ODE> problem,
    AppendState<ODE> append_state,
    Time step)
    : equation(std::move(problem.equation)),
      current_startup_state(*problem.initial_state),
      append_state(std::move(append_state)),
      step(std::move(step)) {
  CHECK_EQ(problem.initial_state->positions.size(),
           problem.initial_state->velocities.size());

  previous_steps.emplace_back();
  FillStepFromSystemState(
      equation, *current_startup_state, previous_steps.back());
}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::StartupSolve(
    Instant const& t_final,
    Instance& instance) const {
  auto const& equation = instance.equation;
  auto const& previous_steps = instance.previous_steps;
  Time const& step = instance.step;
  Time const startup_step = step / startup_step_divisor;

  CHECK(!previous_steps.empty());
  CHECK_LT(previous_steps.size(), order_);

  int startup_step_index = 0;
  auto const startup_append_state =
      [&instance, &startup_step_index](typename ODE::SystemState const& state) {
        CHECK(instance.current_startup_state);
        instance.current_startup_state = state;
        // The startup integrator has a smaller step.  We do not record all the
        // states it computes, but only those that are a multiple of the main
        // integrator step.
        if (++startup_step_index % startup_step_divisor == 0) {
          CHECK_LT(instance.previous_steps.size(), order_);
          instance.previous_steps.emplace_back();
          instance.append_state(state);
          FillStepFromSystemState(instance.equation,
                                  *instance.current_startup_state,
                                  instance.previous_steps.back());
        }
      };

  CHECK(instance.current_startup_state);
  typename ODE::SystemState const& current_startup_state =
      *instance.current_startup_state;
  auto const startup_instance =
      startup_integrator_.NewInstance({equation, &current_startup_state},
                                      startup_append_state,
                                      startup_step);

  startup_integrator_.Solve(
      std::min(current_startup_state.time.value +
                   (order_ - previous_steps.size()) * step + step / 2.0,
               t_final),
      startup_instance.get());

  CHECK_LE(previous_steps.size(), order_);
}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::
FillStepFromSystemState(ODE const& equation,
                        typename ODE::SystemState const& state,
                        Step& step) {
  step.time = state.time;
  for (auto const& position : state.positions) {
    step.displacements.push_back(position.value - Position());
  }
  step.accelerations.resize(step.displacements.size());
  equation.compute_acceleration(step.time.value,
                                step.displacements,
                                &step.accelerations);
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8A() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8A,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -2.0, 2.0, -2.0, 2.0},
      {0.0, 22081.0, -29418.0, 75183.0, -75212.0},
      15120.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8B() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8B,
      BlanesMoan2002SRKN14A<Position>(),
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
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -2.0, 2.0, -1.0, 0.0},
      {0.0, 17671.0, -23622.0, 61449.0, -50516.0},
      12096.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 10> const&
QuinlanTremaine1990Order10() {
  static SymmetricLinearMultistepIntegrator<Position, 10> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_10,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -1.0, 1.0, -1.0, 1.0, -2.0},
      {0.0, 399187.0, -485156.0, 2391436.0, -2816732.0, 4651330.0},
      241920.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 12> const&
QuinlanTremaine1990Order12() {
  static SymmetricLinearMultistepIntegrator<Position, 12> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_12,
      BlanesMoan2002SRKN14A<Position>(),
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
      BlanesMoan2002SRKN14A<Position>(),
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
