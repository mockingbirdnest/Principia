#pragma once

#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "geometry/serialization.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"    // 🧙 For _symplectic_runge_kutta_nyström_integrator.  // NOLINT

namespace principia {
namespace integrators {
namespace _symmetric_linear_multistep_integrator {
namespace internal {

using namespace principia::geometry::_serialization;

int const startup_step_divisor = 16;

template<typename Method, typename ODE_>
absl::Status
SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Solve(
    Instant const& t_final) {
  using Position = typename ODE::DependentVariable;
  using Displacement = typename ODE::DependentVariableDifference;
  using Acceleration = typename ODE::DependentVariableDerivative2;
  using DoubleDisplacement = DoublePrecision<Displacement>;
  using DoubleDisplacements = std::vector<DoubleDisplacement>;
  using DoublePosition = DoublePrecision<Position>;

  auto const& α = integrator_.α_;
  auto const& β_numerator = integrator_.β_numerator_;
  auto const& β_denominator = integrator_.β_denominator_;

  auto& current_state = this->current_state_;
  auto& append_state = this->append_state_;
  auto const& step = this->step_;
  auto const& equation = this->equation_;

  if (!starter_.started()) {
    starter_.Solve(t_final);

    // If `t_final` is not large enough, we may not have generated enough
    // points.  Bail out, we'll continue the next time `Solve` is called.
    if (!starter_.started()) {
      return absl::OkStatus();
    }
  }
  auto const& previous_steps = starter_.previous_steps();

  // Argument checks.
  int const dimension = previous_steps.back().displacements.size();

  // Time step.
  CHECK_LT(Time(), step);
  Time const& h = step;
  // Current time.
  DoublePrecision<Instant> t = previous_steps.back().time;
  // Order.
  int const k = order;

  absl::Status status;
  std::vector<Position> positions(dimension);

  DoubleDisplacements Σⱼ_minus_αⱼ_qⱼ(dimension);
  std::vector<Acceleration> Σⱼ_βⱼ_numerator_aⱼ(dimension);
  while (h <= (t_final - t.value) - t.error) {
    // We take advantage of the symmetry to iterate on the list of previous
    // steps from both ends.
    auto front_it = previous_steps.begin();
    auto back_it = previous_steps.rbegin();

    // This block corresponds to j = 0.  We must not pair it with j = k.
    {
      DoubleDisplacements const& qⱼ = front_it->displacements;
      std::vector<Acceleration> const& aⱼ = front_it->accelerations;
      double const αⱼ = α[0];
      double const βⱼ_numerator = β_numerator[0];
      for (int d = 0; d < dimension; ++d) {
        Σⱼ_minus_αⱼ_qⱼ[d] = Scale(-αⱼ, qⱼ[d]);
        Σⱼ_βⱼ_numerator_aⱼ[d] = βⱼ_numerator * aⱼ[d];
      }
      ++front_it;
    }
    // The generic value of j, paired with k - j.
    for (int j = 1; j < k / 2; ++j) {
      DoubleDisplacements const& qⱼ = front_it->displacements;
      DoubleDisplacements const& qₖ₋ⱼ = back_it->displacements;
      std::vector<Acceleration> const& aⱼ = front_it->accelerations;
      std::vector<Acceleration> const& aₖ₋ⱼ = back_it->accelerations;
      double const αⱼ = α[j];
      double const βⱼ_numerator = β_numerator[j];
      for (int d = 0; d < dimension; ++d) {
        Σⱼ_minus_αⱼ_qⱼ[d] -= Scale(αⱼ, qⱼ[d]);
        Σⱼ_minus_αⱼ_qⱼ[d] -= Scale(αⱼ, qₖ₋ⱼ[d]);
        Σⱼ_βⱼ_numerator_aⱼ[d] += βⱼ_numerator * (aⱼ[d] + aₖ₋ⱼ[d]);
      }
      ++front_it;
      ++back_it;
    }
    // This block corresponds to j = k / 2.  We must not pair it with j = k / 2.
    {
      DoubleDisplacements const& qⱼ = front_it->displacements;
      std::vector<Acceleration> const& aⱼ = front_it->accelerations;
      double const αⱼ = α[k / 2];
      double const βⱼ_numerator = β_numerator[k / 2];
      for (int d = 0; d < dimension; ++d) {
        Σⱼ_minus_αⱼ_qⱼ[d] -= Scale(αⱼ, qⱼ[d]);
        Σⱼ_βⱼ_numerator_aⱼ[d] += βⱼ_numerator * aⱼ[d];
      }
    }

    // Create a new step in the instance.
    t.Increment(h);
    Step current_step{.time = t};
    current_step.accelerations.resize(dimension);
    current_step.displacements.reserve(dimension);

    // Fill the new step.  We skip the division by αₖ as it is equal to 1.0.
    double const αₖ = α[0];
    DCHECK_EQ(αₖ, 1.0);
    for (int d = 0; d < dimension; ++d) {
      DoubleDisplacement& current_displacement = Σⱼ_minus_αⱼ_qⱼ[d];
      current_displacement.Increment(h * h *
                                     Σⱼ_βⱼ_numerator_aⱼ[d] / β_denominator);
      current_step.displacements.push_back(current_displacement);
      DoublePosition const current_position =
          DoublePosition() + current_displacement;
      positions[d] = current_position.value;
      current_state.positions[d] = current_position;
    }
    status.Update(
        equation.compute_acceleration(t.value,
                                      positions,
                                      current_step.accelerations));
    starter_.Push(std::move(current_step));

    ComputeVelocityUsingCohenHubbardOesterwinter();

    // Inform the caller of the new state.
    current_state.time = t;
    append_state(current_state);
    RETURN_IF_STOPPED;  // After the state has been updated.
    if (absl::IsAborted(status)) {
      return status;
    }
  }

  return status;
}

template<typename Method, typename ODE_>
SymmetricLinearMultistepIntegrator<Method, ODE_> const&
SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::integrator()
    const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
void SymmetricLinearMultistepIntegrator<Method, ODE_>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  FixedStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension =
      message
          ->MutableExtension(
              serialization::FixedStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::SymmetricLinearMultistepIntegratorInstance::
                  extension);
  starter_.template WriteToMessage<
      serialization::SymmetricLinearMultistepIntegratorInstance>(extension);
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<
    typename SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance>>
SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::ReadFromMessage(
    serialization::SymmetricLinearMultistepIntegratorInstance const& extension,
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    SymmetricLinearMultistepIntegrator const& integrator)
  requires serializable<typename ODE_::DependentVariable> {
  auto instance = std::unique_ptr<Instance>(new Instance(problem,
                                                         append_state,
                                                         step,
                                                         integrator));
  instance->starter_.FillFromMessage(extension);
  return std::move(instance);
}

template<typename Method, typename ODE_>
void SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Step::
WriteToMessage(
    not_null<serialization::SymmetricLinearMultistepIntegratorInstance::
                  Step*> const message) const {
  using AccelerationSerializer = QuantityOrMultivectorSerializer<
      typename ODE::DependentVariableDerivative2,
      serialization::SymmetricLinearMultistepIntegratorInstance::Step::
          Acceleration>;
  for (auto const& displacement : displacements) {
    displacement.WriteToMessage(message->add_displacements());
  }
  for (auto const& acceleration : accelerations) {
    AccelerationSerializer::WriteToMessage(acceleration,
                                           message->add_accelerations());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename Method, typename ODE_>
typename SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Step
SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Step::
ReadFromMessage(
    serialization::SymmetricLinearMultistepIntegratorInstance::Step const&
        message)
  requires serializable<typename ODE_::DependentVariable> {
  using AccelerationSerializer = QuantityOrMultivectorSerializer<
      typename ODE::DependentVariableDerivative2,
      serialization::SymmetricLinearMultistepIntegratorInstance::Step::
          Acceleration>;
  Step step;
  for (auto const& displacement : message.displacements()) {
    step.displacements.push_back(
        DoublePrecision<typename ODE::DependentVariableDifference>::
            ReadFromMessage(displacement));
  }
  for (auto const& acceleration : message.accelerations()) {
    step.accelerations.push_back(
        AccelerationSerializer::ReadFromMessage(acceleration));
  }
  step.time = DoublePrecision<Instant>::ReadFromMessage(message.time());
  return step;
}

template<typename Method, typename ODE_>
void SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Starter::
FillStepFromState(ODE const& equation,
                  typename ODE::State const& state,
                  Step& step) const {
  using Position = typename ODE::DependentVariable;
  typename ODE::DependentVariables positions;
  step.time = state.time;
  for (auto const& position : state.positions) {
    step.displacements.push_back(position - DoublePrecision<Position>());
    positions.push_back(position.value);
  }
  step.accelerations.resize(step.displacements.size());
  // Ignore the status here.  We are merely computing the acceleration to store
  // it, not to advance an integrator.
  equation.compute_acceleration(step.time.value, positions, step.accelerations)
      .IgnoreError();
}

template<typename Method, typename ODE_>
typename ODE_::IndependentVariable
SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Starter::
independent_variable() const {
  return this->instance().state().time.value;
}

template<typename Method, typename ODE_>
SymmetricLinearMultistepIntegrator<Method, ODE_>::Instance::Instance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    SymmetricLinearMultistepIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      starter_(integrator.startup_integrator_, startup_step_divisor, this),
      integrator_(integrator) {}

template<typename Method, typename ODE_>
void SymmetricLinearMultistepIntegrator<Method, ODE_>::
Instance::ComputeVelocityUsingCohenHubbardOesterwinter() {
  using Displacement = typename ODE::DependentVariableDifference;
  using Velocity = typename ODE::DependentVariableDerivative;
  using Acceleration = typename ODE::DependentVariableDerivative2;

  // For the computation of the velocity we use a formula similar to that of
  // [CHO73].  More specifically, we use the coefficients η from
  // cohen_hubbard_oesterwinter.wl.
  auto const& cohen_hubbard_oesterwinter =
      integrator_.cohen_hubbard_oesterwinter_;

  auto const& previous_steps = starter_.previous_steps();
  int const dimension = previous_steps.back().displacements.size();
  auto& current_state = this->current_state_;
  auto const& step = this->step_;

  current_state.velocities.reserve(dimension);
  for (int d = 0; d < dimension; ++d) {
    DoublePrecision<Velocity>& velocity = current_state.velocities[d];
    auto it = previous_steps.rbegin();

    // Compute the displacement difference using double precision.
    DoublePrecision<Displacement> displacement_change =
        it->displacements[d] - std::next(it)->displacements[d];
    velocity = DoublePrecision<Velocity>(
        (displacement_change.value + displacement_change.error) / step);

    Acceleration weighted_accelerations;
    for (int i = 0;
         i < cohen_hubbard_oesterwinter.numerators.size();
         ++i, ++it) {
      weighted_accelerations +=
          cohen_hubbard_oesterwinter.numerators[i] * it->accelerations[d];
    }

    velocity.value +=
        weighted_accelerations * step / cohen_hubbard_oesterwinter.denominator;
  }
}

template<typename Method, typename ODE_>
SymmetricLinearMultistepIntegrator<Method, ODE_>::
SymmetricLinearMultistepIntegrator(
    FixedStepSizeIntegrator<ODE> const& startup_integrator)
    : startup_integrator_(startup_integrator),
      cohen_hubbard_oesterwinter_(CohenHubbardOesterwinterOrder<order>()) {
  CHECK_EQ(α_[0], 1.0);
  CHECK_EQ(β_numerator_[0], 0.0);
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
SymmetricLinearMultistepIntegrator<Method, ODE_>::NewInstance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step) const {
  // Cannot use `make_not_null_unique` because the constructor of `Instance` is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename ODE_>
void SymmetricLinearMultistepIntegrator<Method, ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal

template<typename Method, typename ODE_>
internal::SymmetricLinearMultistepIntegrator<Method, ODE_> const&
SymmetricLinearMultistepIntegrator() {
  static_assert(
      std::is_base_of<_methods::SymmetricLinearMultistep, Method>::value,
      "Method must be derived from SymmetricLinearMultistep");
  // TODO(phl): Someday, and that day may never come, I will call upon you to
  // expose the startup integrator to the clients.  But until that day, accept
  // this Blanes-Moan integrator as a gift.
  static internal::SymmetricLinearMultistepIntegrator<Method, ODE_> const
      integrator(_symplectic_runge_kutta_nyström_integrator::
                     SymplecticRungeKuttaNyströmIntegrator<
                         _methods::BlanesMoan2002SRKN14A,
                         ODE_>());
  return integrator;
}

}  // namespace _symmetric_linear_multistep_integrator
}  // namespace integrators
}  // namespace principia

#undef PRINCIPIA_USE_COHEN_HUBBARD_OESTERWINTER
