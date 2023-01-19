#pragma once

#include "integrators/explicit_linear_multistep_integrator.hpp"

#include <algorithm>
#include <list>
#include <utility>
#include <vector>

#include "base/jthread.hpp"
#include "geometry/serialization.hpp"
#include "integrators/explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"

namespace principia {
namespace integrators {
namespace internal_explicit_linear_multistep_integrator {

using base::make_not_null_unique;
using geometry::QuantityOrMultivectorSerializer;

int const startup_step_divisor = 16;

template<typename Method, typename ODE_>
absl::Status
ExplicitLinearMultistepIntegrator<Method, ODE_>::Instance::Solve(
    IndependentVariable const& s_final) {
  auto const& α = integrator_.α_;
  auto const& β_numerator = integrator_.β_numerator_;
  auto const& β_denominator = integrator_.β_denominator_;

  auto& current_state = this->current_state_;
  auto& append_state = this->append_state_;
  auto const& step = this->step_;
  auto const& equation = this->equation_;

  if (!starter_.started()) {
    starter_.Solve(s_final);

    // If |s_final| is not large enough, we may not have generated enough
    // points.  Bail out, we'll continue the next time |Solve| is called.
    if (!starter_.started()) {
      return absl::OkStatus();
    }
  }
  auto const& previous_steps = starter_.previous_steps();

  // Independent variable step.
  CHECK_LT(IndependentVariableDifference(), step);
  IndependentVariableDifference const& h = step;
  // Current independent variable.
  DoublePrecision<IndependentVariable> s = previous_steps.back().s;
  // Order.
  int const k = order;

  absl::Status status;

  DoubleDependentVariables Σⱼ_minus_αⱼ_yⱼ;
  DependentVariableDerivatives Σⱼ_βⱼ_numerator_fⱼ;
  while (h <= (s_final - s.value) - s.error) {
    // Scan the previous steps from the most recent to the oldest.  That's how
    // the Adams-Bashforth formulæ are typically written.
    auto it = previous_steps.end();

    for (int j = 1; j <= k; ++j) {
      --it;
      DoubleDependentVariables const& yⱼ = it->y;
      DependentVariableDerivatives const& fⱼ = it->yʹ;
      double const αⱼ = α[j];
      double const βⱼ_numerator = β_numerator[j];
      Σⱼ_minus_αⱼ_yⱼ -= αⱼ * yⱼ;
      Σⱼ_βⱼ_numerator_fⱼ += βⱼ_numerator * fⱼ;
    }

    // Create a new step in the instance.
    s.Increment(h);
    Step current_step{.s = s};

    // Fill the new step.  We skip the division by αₖ as it is equal to 1.0.
    double const αₖ = α[0];
    DCHECK_EQ(αₖ, 1.0);
    Σⱼ_minus_αⱼ_yⱼ.Increment(h * Σⱼ_βⱼ_numerator_fⱼ / β_denominator);
    current_step.y = Σⱼ_minus_αⱼ_yⱼ;
    current_state.y = Σⱼ_minus_αⱼ_yⱼ.value;
    termination_condition::UpdateWithAbort(
        equation.compute_acceleration(s.value,
                                      current_state.y,
                                      current_step.yʹ),
        status);
    starter_.Push(std::move(current_step));

    // Inform the caller of the new state.
    RETURN_IF_STOPPED;
    current_state.s = s;
    append_state(current_state);
    if (absl::IsAborted(status)) {
      return status;
    }
  }

  return status;
}

template<typename Method, typename ODE_>
ExplicitLinearMultistepIntegrator<Method, ODE_> const&
ExplicitLinearMultistepIntegrator<Method, ODE_>::Instance::integrator()
    const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ExplicitLinearMultistepIntegrator<Method, ODE_>::Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
void ExplicitLinearMultistepIntegrator<Method, ODE_>::Instance::Starter::
FillStepFromState(ODE const& equation,
                  typename ODE::State const& state,
                  Step& step) const {
  using Position = typename ODE::DependentVariable;
  typename ODE::DependentVariables positions;
  step.s = state.s;
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
ExplicitLinearMultistepIntegrator<Method, ODE_>::Instance::Instance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    IndependentVariableDifference const& step,
    ExplicitLinearMultistepIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      starter_(integrator.startup_integrator_, startup_step_divisor, this),
      integrator_(integrator) {}

template<typename Method, typename ODE_>
ExplicitLinearMultistepIntegrator<Method, ODE_>::
ExplicitLinearMultistepIntegrator(
    FixedStepSizeIntegrator<ODE> const& startup_integrator)
    : startup_integrator_(startup_integrator) {
  CHECK_EQ(α_[0], 1.0);
  CHECK_EQ(β_numerator_[0], 0.0);
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ExplicitLinearMultistepIntegrator<Method, ODE_>::NewInstance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    IndependentVariableDifference const& step) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename ODE_>
void ExplicitLinearMultistepIntegrator<Method, ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal_explicit_linear_multistep_integrator

template<typename Method, typename ODE_>
internal_explicit_linear_multistep_integrator::
    ExplicitLinearMultistepIntegrator<Method, ODE_> const&
ExplicitLinearMultistepIntegrator() {
  static_assert(
      std::is_base_of<methods::ExplicitLinearMultistep, Method>::value,
      "Method must be derived from ExplicitLinearMultistep");
  // TODO(phl): Someday, and that day may never come, I will call upon you to
  // expose the startup integrator to the clients.  But until that day, accept
  // this Runge-Kutta integrator as a gift.
  static internal_explicit_linear_multistep_integrator::
      ExplicitLinearMultistepIntegrator<Method, ODE_> const integrator(
          ExplicitRungeKuttaIntegrator<methods::Kutta1901Vσ1, ODE_>());
  return integrator;
}

}  // namespace integrators
}  // namespace principia

#undef PRINCIPIA_USE_COHEN_HUBBARD_OESTERWINTER
