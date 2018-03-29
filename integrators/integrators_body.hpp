
#pragma once

#include "integrators/integrators.hpp"

#include <limits>
#include <string>

#include "base/macros.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
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
      current_state_(problem.initial_state),
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
typename ODE_::SystemState const& Integrator<ODE_>::Instance::state() const {
  return current_state_;
}

template<typename ODE_>
void Integrator<ODE_>::Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  current_state_.WriteToMessage(message->mutable_current_state());
}

template<typename ODE_>
Integrator<ODE_>::Instance::Instance() : equation_() {}

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
  IntegrationProblem<ODE> problem;
  problem.equation = equation;
  problem.initial_state =
      ODE::SystemState::ReadFromMessage(message.current_state());

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
      return SymplecticRungeKuttaNyströmIntegrator<
          methods::BlanesMoan2002SRKN6B,
          typename ODE::Position>();
    case FSSI::BLANES_MOAN_2002_SRKN_11B:
      return SymplecticRungeKuttaNyströmIntegrator<
          methods::BlanesMoan2002SRKN11B,
          typename ODE::Position>();
    case FSSI::BLANES_MOAN_2002_SRKN_14A:
      return SymplecticRungeKuttaNyströmIntegrator<
          methods::BlanesMoan2002SRKN14A,
          typename ODE::Position>();
    case FSSI::MCLACHLAN_1995_SB3A_4:
      return SymplecticRungeKuttaNyströmIntegrator<methods::McLachlan1995SB3A4,
                                                   typename ODE::Position>();
    case FSSI::MCLACHLAN_1995_SB3A_5:
      return SymplecticRungeKuttaNyströmIntegrator<methods::McLachlan1995SB3A5,
                                                   typename ODE::Position>();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL:
      return SymplecticRungeKuttaNyströmIntegrator<
          methods::McLachlanAtela1992Order4Optimal,
          typename ODE::Position>();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL:
      return SymplecticRungeKuttaNyströmIntegrator<
          methods::McLachlanAtela1992Order5Optimal,
          typename ODE::Position>();
    case FSSI::OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13:
      return SymplecticRungeKuttaNyströmIntegrator<
          methods::OkunborSkeel1994Order6Method13,
          typename ODE::Position>();
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

template<typename Equation>
FixedStepSizeIntegrator<Equation> const&
ParseFixedStepSizeIntegrator(std::string const& integrator_kind) {
  serialization::FixedStepSizeIntegrator::Kind kind;
  CHECK(serialization::FixedStepSizeIntegrator::Kind_Parse(integrator_kind,
                                                           &kind))
      << "'" << integrator_kind
      << "' is not a valid FixedStepSizeIntegrator.Kind";
  serialization::FixedStepSizeIntegrator message;
  message.set_kind(kind);
  return FixedStepSizeIntegrator<Equation>::ReadFromMessage(message);
}

template<typename ODE_>
AdaptiveStepSizeIntegrator<ODE_>::Parameters::Parameters(
    Time const first_time_step,
    double const safety_factor,
    std::int64_t const max_steps,
    bool const last_step_is_exact)
    : first_time_step(first_time_step),
      safety_factor(safety_factor),
      max_steps(max_steps),
      last_step_is_exact(last_step_is_exact) {}

template<typename ODE_>
AdaptiveStepSizeIntegrator<ODE_>::Parameters::Parameters(
    Time const first_time_step,
    double const safety_factor)
    : Parameters(first_time_step,
                 safety_factor,
                 /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
                 /*last_step_is_exact=*/true) {}

template<typename ODE_>
void AdaptiveStepSizeIntegrator<ODE_>::Parameters::WriteToMessage(
    not_null<serialization::AdaptiveStepSizeIntegratorInstance::
                 Parameters*> const message) const {
  first_time_step.WriteToMessage(message->mutable_first_time_step());
  message->set_safety_factor(safety_factor);
  message->set_max_steps(max_steps);
  message->set_last_step_is_exact(last_step_is_exact);
}

template<typename ODE_>
typename AdaptiveStepSizeIntegrator<ODE_>::Parameters
AdaptiveStepSizeIntegrator<ODE_>::Parameters::ReadFromMessage(
    serialization::AdaptiveStepSizeIntegratorInstance::Parameters const&
        message) {
  bool const is_pre_cartan = !message.has_last_step_is_exact();
  Parameters result(Time::ReadFromMessage(message.first_time_step()),
                    message.safety_factor(),
                    message.max_steps(),
                    is_pre_cartan ? true : message.last_step_is_exact());
  return result;
}

template<typename ODE_>
void AdaptiveStepSizeIntegrator<ODE_>::Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  Integrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension = message->MutableExtension(
      serialization::AdaptiveStepSizeIntegratorInstance::extension);
  parameters_.WriteToMessage(extension->mutable_parameters());
  time_step_.WriteToMessage(extension->mutable_time_step());
  extension->set_first_use(first_use_);
  integrator().WriteToMessage(extension->mutable_integrator());
}

template<typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
AdaptiveStepSizeIntegrator<ODE_>::Instance::ReadFromMessage(
    serialization::IntegratorInstance const& message,
    ODE const& equation,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio) {
  IntegrationProblem<ODE> problem;
  problem.equation = equation;
  problem.initial_state =
      ODE::SystemState::ReadFromMessage(message.current_state());

  CHECK(message.HasExtension(
      serialization::AdaptiveStepSizeIntegratorInstance::extension))
      << "Not an adaptive-step integrator instance " << message.DebugString();
  auto const& extension = message.GetExtension(
      serialization::AdaptiveStepSizeIntegratorInstance::extension);
  auto parameters =
      Parameters::ReadFromMessage(extension.parameters());
  AdaptiveStepSizeIntegrator const& integrator =
      AdaptiveStepSizeIntegrator::ReadFromMessage(extension.integrator());

  // TODO(phl): We would really like this function to return a pointer to an
  // instance from this class, not an instance from |Integrator|.  Unfortunately
  // this confuses Visual Studio 2015...
  auto instance = integrator.ReadFromMessage(
      extension, problem, append_state, tolerance_to_error_ratio, parameters);

  Instance* const down_cast_instance = dynamic_cast<Instance*>(&*instance);
  bool const is_pre_cartan = !extension.has_time_step();
  if (is_pre_cartan) {
    down_cast_instance->time_step_ = parameters.first_time_step;
  } else {
    down_cast_instance->time_step_ =
        Time::ReadFromMessage(extension.time_step());
    down_cast_instance->first_use_ = extension.first_use();
  }

  return instance;
}

template<typename ODE_>
AdaptiveStepSizeIntegrator<ODE_>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters)
    : Integrator<ODE>::Instance(problem, append_state),
      tolerance_to_error_ratio_(tolerance_to_error_ratio),
      parameters_(parameters),
      time_step_(parameters.first_time_step) {
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

template<typename Equation>
AdaptiveStepSizeIntegrator<Equation> const& ParseAdaptiveStepSizeIntegrator(
    std::string const& integrator_kind) {
  serialization::AdaptiveStepSizeIntegrator::Kind kind;
  CHECK(serialization::AdaptiveStepSizeIntegrator::Kind_Parse(integrator_kind,
                                                              &kind))
      << "'" << integrator_kind
      << "' is not a valid AdaptiveStepSizeIntegrator.Kind";
  serialization::AdaptiveStepSizeIntegrator message;
  message.set_kind(kind);
  return AdaptiveStepSizeIntegrator<Equation>::ReadFromMessage(message);
}

}  // namespace internal_integrators
}  // namespace integrators
}  // namespace principia
