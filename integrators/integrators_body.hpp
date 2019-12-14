
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

// A case branch in a switch on the serialized integrator |kind|.  It determines
// the |method| type from the |kind| and calls |action| on it.  |action| must be
// a 1-argument macro.
#define PRINCIPIA_INTEGRATOR_CASE(kind, method, action) \
  case serialization::FixedStepSizeIntegrator::kind: {  \
    action(method);                                     \
    break;                                              \
  }

// All the case branches in a switch on the serialized integrator kind.
// Depending on the nature of the integrator, each case branch calls one of the
// given actions, which must be 1-argument macros.
#define PRINCIPIA_INTEGRATOR_CASES(slms_action, sprk_action, srkn_action)      \
  PRINCIPIA_INTEGRATOR_CASE(BLANES_MOAN_2002_S6,                               \
                            BlanesMoan2002S6,                                  \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(BLANES_MOAN_2002_S10,                              \
                            BlanesMoan2002S10,                                 \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(BLANES_MOAN_2002_SRKN_6B,                          \
                            BlanesMoan2002SRKN6B,                              \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(BLANES_MOAN_2002_SRKN_11B,                         \
                            BlanesMoan2002SRKN11B,                             \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(BLANES_MOAN_2002_SRKN_14A,                         \
                            BlanesMoan2002SRKN14A,                             \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(CANDY_ROZMUS_1991_FOREST_RUTH_1990,                \
                            CandyRozmus1991ForestRuth1990,                     \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_S2,                                 \
                            McLachlan1995S2,                                   \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_S4,                                 \
                            McLachlan1995S4,                                   \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_S5,                                 \
                            McLachlan1995S5,                                   \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_SB3A_4,                             \
                            McLachlan1995SB3A4,                                \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_SB3A_5,                             \
                            McLachlan1995SB3A5,                                \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_SS5,                                \
                            McLachlan1995SS5,                                  \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_SS9,                                \
                            McLachlan1995SS9,                                  \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_SS15,                               \
                            McLachlan1995SS15,                                 \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_1995_SS17,                               \
                            McLachlan1995SS17,                                 \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_ATELA_1992_ORDER_2_OPTIMAL,              \
                            McLachlanAtela1992Order2Optimal,                   \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_ATELA_1992_ORDER_3_OPTIMAL,              \
                            McLachlanAtela1992Order3Optimal,                   \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL,              \
                            McLachlanAtela1992Order4Optimal,                   \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL,              \
                            McLachlanAtela1992Order5Optimal,                   \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(NEWTON_DELAMBRE_STORMER_VERLET_LEAPFROG,           \
                            NewtonDelambreStørmerVerletLeapfrog,               \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13,              \
                            OkunborSkeel1994Order6Method13,                    \
                            srkn_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(QUINLAN_1999_ORDER_8A,                             \
                            Quinlan1999Order8A,                                \
                            slms_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(QUINLAN_1999_ORDER_8B,                             \
                            Quinlan1999Order8B,                                \
                            slms_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(QUINLAN_TREMAINE_1990_ORDER_8,                     \
                            QuinlanTremaine1990Order8,                         \
                            slms_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(QUINLAN_TREMAINE_1990_ORDER_10,                    \
                            QuinlanTremaine1990Order10,                        \
                            slms_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(QUINLAN_TREMAINE_1990_ORDER_12,                    \
                            QuinlanTremaine1990Order12,                        \
                            slms_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(QUINLAN_TREMAINE_1990_ORDER_14,                    \
                            QuinlanTremaine1990Order14,                        \
                            slms_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(RUTH_1983,                                         \
                            Ruth1983,                                          \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(SUZUKI_1990,                                       \
                            鈴木1990,                                           \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_6A,                             \
                            吉田1990Order6A,                                    \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_6B,                             \
                            吉田1990Order6B,                                    \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_6C,                             \
                            吉田1990Order6C,                                    \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_8A,                             \
                            吉田1990Order8A,                                    \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_8B,                             \
                            吉田1990Order8B,                                    \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_8C,                             \
                            吉田1990Order8C,                                    \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_8D,                             \
                            吉田1990Order8D,                                    \
                            sprk_action)                                       \
  PRINCIPIA_INTEGRATOR_CASE(YOSHIDA_1990_ORDER_8E,                             \
                            吉田1990Order8E,                                    \
                            sprk_action)

namespace principia {
namespace integrators {
namespace internal_integrators {

template<typename Integrator>
not_null<std::unique_ptr<typename Integrator::Instance>>
ReadSlmsInstanceFromMessage(
    serialization::FixedStepSizeIntegratorInstance const& message,
    IntegrationProblem<typename Integrator::ODE> const& problem,
    typename Integrator::AppendState const& append_state,
    Time const& step,
    Integrator const& integrator) {
  CHECK(message.HasExtension(
      serialization::SymmetricLinearMultistepIntegratorInstance::extension))
      << message.DebugString();
  auto const& extension = message.GetExtension(
      serialization::SymmetricLinearMultistepIntegratorInstance::extension);
  return typename Integrator::Instance::ReadFromMessage(extension,
                                                        problem,
                                                        append_state,
                                                        step,
                                                        integrator);
}

template<typename Integrator>
not_null<std::unique_ptr<typename Integrator::Instance>>
ReadSrknInstanceFromMessage(
    serialization::FixedStepSizeIntegratorInstance const& message,
    IntegrationProblem<typename Integrator::ODE> const& problem,
    typename Integrator::AppendState const& append_state,
    Time const& step,
    Integrator const& integrator) {
  CHECK(message.HasExtension(
      serialization::SymplecticRungeKuttaNystromIntegratorInstance::extension))
      << message.DebugString();
  auto const& extension = message.GetExtension(
      serialization::SymplecticRungeKuttaNystromIntegratorInstance::extension);
  return typename Integrator::Instance::ReadFromMessage(extension,
                                                        problem,
                                                        append_state,
                                                        step,
                                                        integrator);
}

// We do not deserialize an SPRK per se, but only when it is converted to an
// SRKN.  The reason is that an SPRK is for a different kind of equation than
// an SRKN, so the two would return different types.  We avoid replicating code
// at the expense of some template complexity.
template<typename Integrator, typename Result>
struct SprkAsSrknConstructor;

template<typename Integrator>
struct SprkAsSrknConstructor<Integrator, Integrator const&> {
  static Integrator const& Make(
      serialization::FixedStepSizeIntegratorInstance const& message,
      Integrator const& integrator) {
    return integrator;
  }
};

template<typename Integrator, typename Instance>
struct SprkAsSrknConstructor<Integrator, not_null<std::unique_ptr<Instance>>> {
  static not_null<std::unique_ptr<Instance>> Make(
      serialization::FixedStepSizeIntegratorInstance const& message,
      IntegrationProblem<typename Integrator::ODE> const& problem,
      typename Integrator::AppendState const& append_state,
      Time const& step,
      Integrator const& integrator) {
    return ReadSrknInstanceFromMessage(
        message, problem, append_state, step, integrator);
  }
};

template<typename Result,
         typename Position,
         typename Method,
         bool first_same_as_last,
         typename... Args>
Result ReadSprkFromMessage(
    serialization::FixedStepSizeIntegratorInstance const& message,
    Args... args) {
  CHECK(message.integrator().has_composition_method());
  if constexpr (first_same_as_last) {
    switch (message.integrator().composition_method()) {
      case serialization::FixedStepSizeIntegrator::ABA: {
        auto const& integrator = SymplecticRungeKuttaNyströmIntegrator<
            Method,
            serialization::FixedStepSizeIntegrator::ABA,
            Position>();
        using Integrator =
            std::remove_reference_t<std::remove_const_t<decltype(integrator)>>;
        return SprkAsSrknConstructor<Integrator, Result>::Make(
            message, args..., integrator);
      }
      case serialization::FixedStepSizeIntegrator::BAB: {
        auto const& integrator = SymplecticRungeKuttaNyströmIntegrator<
            Method,
            serialization::FixedStepSizeIntegrator::BAB,
            Position>();
        using Integrator =
            std::remove_reference_t<std::remove_const_t<decltype(integrator)>>;
        return SprkAsSrknConstructor<Integrator, Result>::Make(
            message, args..., integrator);
      }
      case serialization::FixedStepSizeIntegrator::BA:
      default:
        LOG(FATAL) << message.DebugString();
        base::noreturn();
    }
  } else {
    CHECK_EQ(serialization::FixedStepSizeIntegrator::BA,
             message.integrator().composition_method());
    auto const& integrator = SymplecticRungeKuttaNyströmIntegrator<
        Method,
        serialization::FixedStepSizeIntegrator::BA,
        Position>();
    using Integrator =
        std::remove_reference_t<std::remove_const_t<decltype(integrator)>>;
    return SprkAsSrknConstructor<Integrator, Result>::Make(
        message, args..., integrator);
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

#define PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SLMS(method)    \
  auto const& integrator =                                            \
      SymmetricLinearMultistepIntegrator<methods::method,             \
                                         typename ODE::Position>();   \
  return ReadSlmsInstanceFromMessage(                                 \
      extension, problem, append_state, step, integrator)

#define PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SPRK(method)   \
  return ReadSprkFromMessage<                                        \
      not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>, \
      typename ODE::Position,                                        \
      methods::method,                                               \
      methods::method::first_same_as_last>(                          \
      extension, problem, append_state, step)

#define PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SRKN(method)     \
  auto const& integrator =                                             \
      SymplecticRungeKuttaNyströmIntegrator<methods::method,           \
                                            typename ODE::Position>(); \
  return ReadSrknInstanceFromMessage(                                  \
      extension, problem, append_state, step, integrator)

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

  switch (extension.integrator().kind()) {
    PRINCIPIA_INTEGRATOR_CASES(
        PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SLMS,
        PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SPRK,
        PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SRKN)
    default:
      LOG(FATAL) << message.DebugString();
      base::noreturn();
  }
}

#undef PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SLMS
#undef PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SPRK
#undef PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_INSTANCE_SPRK

template<typename ODE_>
FixedStepSizeIntegrator<ODE_>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step)
    : Integrator<ODE>::Instance(problem, std::move(append_state)),
      step_(step) {
  CHECK_NE(Time(), step_);
}

#define PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SLMS(method)    \
  return SymmetricLinearMultistepIntegrator<methods::method, \
                                            typename ODE::Position>()

#define PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SPRK(method)  \
  return ReadSprkFromMessage<FixedStepSizeIntegrator<ODE>, \
                             typename ODE::Position,       \
                             methods::method,              \
                             methods::method::first_same_as_last>(message)

#define PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SRKN(method)       \
  return SymplecticRungeKuttaNyströmIntegrator<methods::method, \
                                               typename ODE::Position>()

template<typename ODE_>
FixedStepSizeIntegrator<ODE_> const&
FixedStepSizeIntegrator<ODE_>::ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message) {
  switch (message.kind()) {
    PRINCIPIA_INTEGRATOR_CASES(PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SLMS,
                               PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SPRK,
                               PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SRKN)
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

#undef PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SLMS
#undef PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SPRK
#undef PRINCIPIA_READ_FIXED_STEP_INTEGRATOR_SPRK

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
      template ReadFromMessage<AdaptiveStepSizeIntegrator<ODE_>>(message);
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
#undef PRINCIPIA_CASES
