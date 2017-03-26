
#pragma once

#include "integrators/integrators.hpp"

#include "base/macros.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {
namespace internal_integrators {

template<typename ODE_>
Integrator<ODE_>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state)
    : equation_(problem.equation),
      current_state_(*problem.initial_state),
      append_state_(std::move(append_state)) {
  CHECK_EQ(current_state_.positions.size(),
           current_state_.velocities.size());
}

template<typename ODE_>
DoublePrecision<Instant> const&
Integrator<ODE_>::Instance::time() const {
  return current_state_.time;
}

template<typename ODE_>
void Integrator<ODE_>::Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  current_state_.WriteToMessage(message->mutable_current_state());
}

template<typename ODE_>
Integrator<ODE_>::Instance::Instance() {}

template<typename ODE_>
void FixedStepSizeIntegrator<ODE_>::Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  Integrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension = message->MutableExtension(
      serialization::FixedStepSizeIntegratorInstance::extension);
  step_.WriteToMessage(extension->mutable_step());
  integrator().WriteToMessage(extension->mutable_integrator());
}

template<typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
FixedStepSizeIntegrator<ODE_>::Instance::ReadFromMessage(
    serialization::IntegratorInstance const& message,
    ODE const& equation,
    AppendState const& append_state) {
  auto const current_state =
      ODE::SystemState::ReadFromMessage(message.current_state());
  IntegrationProblem<ODE> problem;
  problem.equation = equation;
  problem.initial_state = &current_state;

  CHECK(message.HasExtension(
      serialization::FixedStepSizeIntegratorInstance::extension))
      << "Not a fixed-step integrator instance " << message.DebugString();
  auto const& extension = message.GetExtension(
      serialization::FixedStepSizeIntegratorInstance::extension);
  Time const step = Time::ReadFromMessage(extension.step());
  FixedStepSizeIntegrator const& integrator =
      FixedStepSizeIntegrator::ReadFromMessage(extension.integrator());

  return integrator.ReadFromMessage(extension, problem, append_state, step);
}

template<typename ODE_>
FixedStepSizeIntegrator<ODE_>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step)
    : Integrator<ODE>::Instance(problem, std::move(append_state)),
      step_(step) {
  CHECK_NE(Time(), step_);
}

template<typename ODE_>
void FixedStepSizeIntegrator<ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> const message) const {
  message->set_kind(kind_);
}

template<typename ODE_>
FixedStepSizeIntegrator<ODE_> const&
FixedStepSizeIntegrator<ODE_>::ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message) {
  using FSSI = serialization::FixedStepSizeIntegrator;
  switch (message.kind()) {
    case FSSI::BLANES_MOAN_2002_SRKN_6B:
      return BlanesMoan2002SRKN6B<typename ODE::Position>();
    case FSSI::BLANES_MOAN_2002_SRKN_11B:
      return BlanesMoan2002SRKN11B<typename ODE::Position>();
    case FSSI::BLANES_MOAN_2002_SRKN_14A:
      return BlanesMoan2002SRKN14A<typename ODE::Position>();
    case FSSI::MCLACHLAN_1995_SB3A_4:
      return McLachlan1995SB3A4<typename ODE::Position>();
    case FSSI::MCLACHLAN_1995_SB3A_5:
      return McLachlan1995SB3A5<typename ODE::Position>();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL:
      return McLachlanAtela1992Order4Optimal<typename ODE::Position>();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL:
      return McLachlanAtela1992Order5Optimal<typename ODE::Position>();
    case FSSI::OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13:
      return OkunborSkeel1994Order6Method13<typename ODE::Position>();
    case FSSI::QUINLAN_1999_ORDER_8A:
      return Quinlan1999Order8A<typename ODE::Position>();
    case FSSI::QUINLAN_1999_ORDER_8B:
      return Quinlan1999Order8B<typename ODE::Position>();
    case FSSI::QUINLAN_TREMAINE_1990_ORDER_8:
      return QuinlanTremaine1990Order8<typename ODE::Position>();
    case FSSI::QUINLAN_TREMAINE_1990_ORDER_10:
      return QuinlanTremaine1990Order10<typename ODE::Position>();
    case FSSI::QUINLAN_TREMAINE_1990_ORDER_12:
      return QuinlanTremaine1990Order12<typename ODE::Position>();
    case FSSI::QUINLAN_TREMAINE_1990_ORDER_14:
      return QuinlanTremaine1990Order14<typename ODE::Position>();
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

template<typename ODE_>
FixedStepSizeIntegrator<ODE_>::FixedStepSizeIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind) : kind_(kind) {}

template<typename ODE_>
void AdaptiveStepSizeIntegrator<ODE_>::Parameters::WriteToMessage(
    not_null<serialization::AdaptiveStepSizeIntegratorInstance::
                 Parameters*> const message) const {
  first_time_step.WriteToMessage(message->mutable_first_time_step());
  message->set_safety_factor(safety_factor);
  message->set_max_steps(max_steps);
}

template<typename ODE_>
typename AdaptiveStepSizeIntegrator<ODE_>::Parameters
AdaptiveStepSizeIntegrator<ODE_>::Parameters::ReadFromMessage(
    serialization::AdaptiveStepSizeIntegratorInstance::Parameters const&
        message) {
  Parameters result;
  result.first_time_step = Time::ReadFromMessage(message.first_time_step());
  result.safety_factor = message.safety_factor();
  result.max_steps = message.max_steps();
  return result;
}

template<typename ODE_>
void AdaptiveStepSizeIntegrator<ODE_>::Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  Integrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension = message->MutableExtension(
      serialization::AdaptiveStepSizeIntegratorInstance::extension);
  parameters_.WriteToMessage(extension->mutable_parameters());
  integrator().WriteToMessage(extension->mutable_integrator());
}

template<typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
AdaptiveStepSizeIntegrator<ODE_>::Instance::ReadFromMessage(
    serialization::IntegratorInstance const& message,
    ODE const& equation,
    AppendState const& append_state,
    typename Parameters::ToleranceToErrorRatio const&
        tolerance_to_error_ratio) {
  auto const current_state =
      ODE::SystemState::ReadFromMessage(message.current_state());
  IntegrationProblem<ODE> problem;
  problem.equation = equation;
  problem.initial_state = &current_state;

  CHECK(message.HasExtension(
      serialization::AdaptiveStepSizeIntegratorInstance::extension))
      << "Not an adaptive-step integrator instance " << message.DebugString();
  auto const& extension = message.GetExtension(
      serialization::AdaptiveStepSizeIntegratorInstance::extension);
  auto parameters =
      Parameters::ReadFromMessage(extension.parameters());
  parameters.tolerance_to_error_ratio = tolerance_to_error_ratio;
  AdaptiveStepSizeIntegrator const& integrator =
      AdaptiveStepSizeIntegrator::ReadFromMessage(extension.integrator());

  return integrator.ReadFromMessage(
      extension, problem, append_state, parameters);
}

template<typename ODE_>
AdaptiveStepSizeIntegrator<ODE_>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Parameters const& parameters)
    : Integrator<ODE>::Instance(problem, std::move(append_state)),
      parameters_(parameters) {
  CHECK_NE(Time(), parameters.first_time_step);
  CHECK_GT(parameters.safety_factor, 0);
  CHECK_LT(parameters.safety_factor, 1);
}

template<typename ODE_>
void AdaptiveStepSizeIntegrator<ODE_>::WriteToMessage(
    not_null<serialization::AdaptiveStepSizeIntegrator*> const message) const {
  message->set_kind(kind_);
}

template<typename ODE_>
AdaptiveStepSizeIntegrator<ODE_> const&
AdaptiveStepSizeIntegrator<ODE_>::ReadFromMessage(
    serialization::AdaptiveStepSizeIntegrator const& message) {
  using ASSI = serialization::AdaptiveStepSizeIntegrator;
  switch (message.kind()) {
    case ASSI::DORMAND_ELMIKKAWY_PRINCE_1986_RKN_434FM:
      return DormandElMikkawyPrince1986RKN434FM<typename ODE::Position>();
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

template<typename ODE_>
AdaptiveStepSizeIntegrator<ODE_>::AdaptiveStepSizeIntegrator(
    serialization::AdaptiveStepSizeIntegrator::Kind const kind) : kind_(kind) {}

}  // namespace internal_integrators
}  // namespace integrators
}  // namespace principia
