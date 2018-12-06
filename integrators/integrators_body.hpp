
#pragma once

#include "integrators/integrators.hpp"

#include <limits>
#include <vector>
#include <string>

#include "base/macros.hpp"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#define PRINCIPIA_CASE_SLMS(kind, method)                      \
  case serialization::FixedStepSizeIntegrator::kind:           \
    return SymmetricLinearMultistepIntegrator<methods::method, \
                                              typename ODE::Position>()

// We do not deserialize an SPRK per se, but only when it is converted to an
// SRKN.  The reason is that an SPRK is for a different kind of equation than
// an SRKN, so the two would return different types.  If we ever need to do
// this we will need to specialize ReadFromMessage somehow.
#define PRINCIPIA_CASE_SPRK(kind, method)                                 \
  case serialization::FixedStepSizeIntegrator::kind:                      \
    do {                                                                  \
      CHECK(message.has_composition_method());                            \
      return SprkAsSrknDeserializer<                                      \
          ODE,                                                            \
          methods::method,                                                \
          methods::method::first_same_as_last>::ReadFromMessage(message); \
    } while (false)

#define PRINCIPIA_CASE_SRKN(kind, method)                         \
  case serialization::FixedStepSizeIntegrator::kind:              \
    return SymplecticRungeKuttaNyströmIntegrator<methods::method, \
                                                 typename ODE::Position>()

namespace principia {
namespace integrators {
namespace internal_integrators {

template<typename ODE, typename Method, bool first_same_as_last>
struct SprkAsSrknDeserializer;

template<typename ODE, typename Method>
struct SprkAsSrknDeserializer<ODE, Method, false> {
  static FixedStepSizeIntegrator<ODE> const& ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message);
};

template<typename ODE, typename Method>
struct SprkAsSrknDeserializer<ODE, Method, true> {
  static FixedStepSizeIntegrator<ODE> const& ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message);
};

template<typename ODE, typename Method>
FixedStepSizeIntegrator<ODE> const&
SprkAsSrknDeserializer<ODE, Method, false>::ReadFromMessage(
    serialization::FixedStepSizeIntegrator const& message) {
  CHECK_EQ(serialization::FixedStepSizeIntegrator::BA,
           message.composition_method());
  return SymplecticRungeKuttaNyströmIntegrator<
      Method,
      serialization::FixedStepSizeIntegrator::BA,
      typename ODE::Position>();
}

template<typename ODE, typename Method>
FixedStepSizeIntegrator<ODE> const&
SprkAsSrknDeserializer<ODE, Method, true>::ReadFromMessage(
    serialization::FixedStepSizeIntegrator const& message) {
  switch (message.composition_method()) {
    case serialization::FixedStepSizeIntegrator::ABA:
      return SymplecticRungeKuttaNyströmIntegrator<
          Method,
          serialization::FixedStepSizeIntegrator::ABA,
          typename ODE::Position>();
    case serialization::FixedStepSizeIntegrator::BAB:
      return SymplecticRungeKuttaNyströmIntegrator<
          Method,
          serialization::FixedStepSizeIntegrator::BAB,
          typename ODE::Position>();
    case serialization::FixedStepSizeIntegrator::BA:
    default:
      LOG(FATAL) << message.DebugString();
      base::noreturn();
  }
}

// The parameters of the |RightHandSideComputation| determines the type of the
// integrator when the method may be used for several kinds of integrators
// (this may not hold in general as we add more integrators).
template<typename RightHandSideComputation>
struct AdaptiveStepSizeIntegratorDeserializer;

template<typename Position, typename Acceleration>
struct AdaptiveStepSizeIntegratorDeserializer<
    std::function<Status(Instant const& t,
                         std::vector<Position> const& positions,
                         std::vector<Acceleration>& accelerations)>> {
  template<typename Integrator>
  static Integrator const& ReadFromMessage(
      serialization::AdaptiveStepSizeIntegrator const& message);
};

template<typename Position, typename Velocity, typename Acceleration>
struct AdaptiveStepSizeIntegratorDeserializer<
    std::function<Status(Instant const& t,
                         std::vector<Position> const& positions,
                         std::vector<Velocity> const& velocities,
                         std::vector<Acceleration>& accelerations)>> {
  template<typename Integrator>
  static Integrator const& ReadFromMessage(
      serialization::AdaptiveStepSizeIntegrator const& message);
};

template<typename Position, typename Acceleration>
template<typename Integrator>
Integrator const& AdaptiveStepSizeIntegratorDeserializer<
    std::function<Status(Instant const& t,
                         std::vector<Position> const& positions,
                         std::vector<Acceleration>& accelerations)>>::
ReadFromMessage(serialization::AdaptiveStepSizeIntegrator const& message) {
  using ASSI = serialization::AdaptiveStepSizeIntegrator;
  switch (message.kind()) {
    case ASSI::DORMAND_ELMIKKAWY_PRINCE_1986_RKN_434FM:
      return EmbeddedExplicitRungeKuttaNyströmIntegrator<
          methods::DormandالمكاوىPrince1986RKN434FM, Position>();
    case ASSI::FINE_1987_RKNG_34:
      return EmbeddedExplicitRungeKuttaNyströmIntegrator<
          methods::Fine1987RKNG34, Position>();
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

template<typename Position, typename Velocity, typename Acceleration>
template<typename Integrator>
Integrator const& AdaptiveStepSizeIntegratorDeserializer<
    std::function<Status(Instant const& t,
                         std::vector<Position> const& positions,
                         std::vector<Velocity> const& velocities,
                         std::vector<Acceleration>& accelerations)>>::
ReadFromMessage(serialization::AdaptiveStepSizeIntegrator const& message) {
  using ASSI = serialization::AdaptiveStepSizeIntegrator;
  switch (message.kind()) {
    case ASSI::FINE_1987_RKNG_34:
      return EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
          methods::Fine1987RKNG34, Position>();
    case ASSI::DORMAND_ELMIKKAWY_PRINCE_1986_RKN_434FM:
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

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
FixedStepSizeIntegrator<ODE_> const&
FixedStepSizeIntegrator<ODE_>::ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message) {
  switch (message.kind()) {
    PRINCIPIA_CASE_SPRK(BLANES_MOAN_2002_S6,
                        BlanesMoan2002S6);
    PRINCIPIA_CASE_SPRK(BLANES_MOAN_2002_S10,
                        BlanesMoan2002S10);
    PRINCIPIA_CASE_SRKN(BLANES_MOAN_2002_SRKN_6B,
                        BlanesMoan2002SRKN6B);
    PRINCIPIA_CASE_SRKN(BLANES_MOAN_2002_SRKN_11B,
                        BlanesMoan2002SRKN11B);
    PRINCIPIA_CASE_SRKN(BLANES_MOAN_2002_SRKN_14A,
                        BlanesMoan2002SRKN14A);
    PRINCIPIA_CASE_SPRK(CANDY_ROZMUS_1991_FOREST_RUTH_1990,
                        CandyRozmus1991ForestRuth1990);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_1995_S2,
                        McLachlan1995S2);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_1995_S4,
                        McLachlan1995S4);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_1995_S5,
                        McLachlan1995S5);
    PRINCIPIA_CASE_SRKN(MCLACHLAN_1995_SB3A_4,
                        McLachlan1995SB3A4);
    PRINCIPIA_CASE_SRKN(MCLACHLAN_1995_SB3A_5,
                        McLachlan1995SB3A5);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_1995_SS5,
                        McLachlan1995SS5);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_1995_SS9,
                        McLachlan1995SS9);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_1995_SS15,
                        McLachlan1995SS15);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_1995_SS17,
                        McLachlan1995SS17);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_ATELA_1992_ORDER_2_OPTIMAL,
                        McLachlanAtela1992Order2Optimal);
    PRINCIPIA_CASE_SPRK(MCLACHLAN_ATELA_1992_ORDER_3_OPTIMAL,
                        McLachlanAtela1992Order3Optimal);
    PRINCIPIA_CASE_SRKN(MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL,
                        McLachlanAtela1992Order4Optimal);
    PRINCIPIA_CASE_SRKN(MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL,
                        McLachlanAtela1992Order5Optimal);
    PRINCIPIA_CASE_SPRK(NEWTON_DELAMBRE_STORMER_VERLET_LEAPFROG,
                        NewtonDelambreStørmerVerletLeapfrog);
    PRINCIPIA_CASE_SRKN(OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13,
                        OkunborSkeel1994Order6Method13);
    PRINCIPIA_CASE_SLMS(QUINLAN_1999_ORDER_8A,
                        Quinlan1999Order8A);
    PRINCIPIA_CASE_SLMS(QUINLAN_1999_ORDER_8B,
                        Quinlan1999Order8B);
    PRINCIPIA_CASE_SLMS(QUINLAN_TREMAINE_1990_ORDER_8,
                        QuinlanTremaine1990Order8);
    PRINCIPIA_CASE_SLMS(QUINLAN_TREMAINE_1990_ORDER_10,
                        QuinlanTremaine1990Order10);
    PRINCIPIA_CASE_SLMS(QUINLAN_TREMAINE_1990_ORDER_12,
                        QuinlanTremaine1990Order12);
    PRINCIPIA_CASE_SLMS(QUINLAN_TREMAINE_1990_ORDER_14,
                        QuinlanTremaine1990Order14);
    PRINCIPIA_CASE_SPRK(RUTH_1983,
                        Ruth1983);
    PRINCIPIA_CASE_SPRK(SUZUKI_1990,
                        鈴木1990);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_6A,
                        吉田1990Order6A);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_6B,
                        吉田1990Order6B);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_6C,
                        吉田1990Order6C);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_8A,
                        吉田1990Order8A);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_8B,
                        吉田1990Order8B);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_8C,
                        吉田1990Order8C);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_8D,
                        吉田1990Order8D);
    PRINCIPIA_CASE_SPRK(YOSHIDA_1990_ORDER_8E,
                        吉田1990Order8E);
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

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
AdaptiveStepSizeIntegrator<ODE_> const&
AdaptiveStepSizeIntegrator<ODE_>::ReadFromMessage(
    serialization::AdaptiveStepSizeIntegrator const& message) {
  return AdaptiveStepSizeIntegratorDeserializer<
      typename ODE_::RightHandSideComputation>::
      ReadFromMessage<AdaptiveStepSizeIntegrator<ODE_>>(message);
}

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

#undef PRINCIPIA_CASE_SLMS
#undef PRINCIPIA_CASE_SPRK
#undef PRINCIPIA_CASE_SRKN
