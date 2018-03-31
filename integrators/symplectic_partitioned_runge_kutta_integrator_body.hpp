
#pragma once

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include "base/mod.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_partitioned_runge_kutta_integrator {

using base::mod;

template<typename Method, typename Position>
Status SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
Instance::Solve(Instant const& t_final) {
  LOG(FATAL) << "I'm sorry Dave.  I'm afraid I can't do that.";
}

template<typename Method, typename Position>
SymplecticPartitionedRungeKuttaIntegrator<Method, Position> const&
SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    DecomposableFirstOrderDifferentialEquation<Position>>::Instance>>
SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename Position>
void SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  FixedStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension =
      message
          ->MutableExtension(
              serialization::FixedStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::SymplecticPartitionedRungeKuttaIntegratorInstance::
                  extension);
}

template<typename Method, typename Position>
SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
SymplecticPartitionedRungeKuttaIntegrator() {
  // TODO(phl): This might be turned into a static_assert.
  if (first_same_as_last) {
    CHECK_EQ(0.0, a_[stages_ - 1]);
  }
  if (time_reversible) {
    CHECK(first_same_as_last);
    for (int i = 0; i < stages_ - 1; ++i) {
      CHECK_EQ(a_[i], a_[stages_ - 2 - i]);
    }
    for (int i = 0; i < stages_; ++i) {
      CHECK_EQ(b_[i], b_[stages_ - 1 - i]);
    }
  }
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    DecomposableFirstOrderDifferentialEquation<Position>>::Instance>>
SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
NewInstance(IntegrationProblem<ODE> const& problem,
            AppendState const& append_state,
            Time const& step) const {
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename Position>
void SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
    WriteToMessage(
        not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    DecomposableFirstOrderDifferentialEquation<Position>>::Instance>>
SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
ReadFromMessage(serialization::FixedStepSizeIntegratorInstance const& message,
                IntegrationProblem<ODE> const& problem,
                AppendState const& append_state,
                Time const& step) const {
  CHECK(message.HasExtension(
      serialization::SymplecticRungeKuttaNystromIntegratorInstance::extension))
      << message.DebugString();

  return std::unique_ptr<typename Integrator<ODE>::Instance>(
      new Instance(problem, append_state, step, *this));
}

}  // namespace internal_symplectic_partitioned_runge_kutta_integrator

template<typename Method, typename Position>
internal_symplectic_partitioned_runge_kutta_integrator::
    SymplecticPartitionedRungeKuttaIntegrator<Method, Position> const&
    SymplecticPartitionedRungeKuttaIntegrator() {
  static_assert(
      std::is_base_of<methods::SymplecticPartitionedRungeKutta, Method>::value,
      "Method must be derived from SymplecticPartitionedRungeKutta");
  static internal_symplectic_partitioned_runge_kutta_integrator::
      SymplecticPartitionedRungeKuttaIntegrator<Method, Position> const
          integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
