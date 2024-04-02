#pragma once

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <memory>

#include "base/mod.hpp"

namespace principia {
namespace integrators {
namespace _symplectic_partitioned_runge_kutta_integrator {
namespace internal {

using namespace principia::base::_mod;

template<typename Method, typename ODE_>
absl::Status SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_>::
Instance::Solve(Instant const& t_final) {
  LOG(FATAL) << "I'm sorry Dave.  I'm afraid I can't do that.";
}

template<typename Method, typename ODE_>
SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_> const&
SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_>::
Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
void SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  LOG(FATAL) << "WriteToMessage NYI";
}

template<typename Method, typename ODE_>
SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_>::
SymplecticPartitionedRungeKuttaIntegrator() {
  static_assert(!first_same_as_last || a_[stages_ - 1] == 0.0);
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

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_>::
NewInstance(InitialValueProblem<ODE> const& problem,
            AppendState const& append_state,
            Time const& step) const {
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename ODE_>
void SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_>::
    WriteToMessage(
        not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal

template<typename Method, typename ODE_>
internal::SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_> const&
SymplecticPartitionedRungeKuttaIntegrator() {
  static_assert(
      std::is_base_of<_methods::SymplecticPartitionedRungeKutta, Method>::value,
      "Method must be derived from SymplecticPartitionedRungeKutta");
  static internal::SymplecticPartitionedRungeKuttaIntegrator<Method, ODE_> const
      integrator;
  return integrator;
}

}  // namespace _symplectic_partitioned_runge_kutta_integrator
}  // namespace integrators
}  // namespace principia
