#pragma once

#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "integrators/ordinary_differential_equations.hpp"

namespace principia {
namespace integrators {
namespace internal_ordinary_differential_equations {

using base::make_not_null_unique;

template <typename DifferentialEquation>
class MockFixedStepSizeIntegrator
    : public FixedStepSizeIntegrator<DifferentialEquation> {
 public:
  using ODE = DifferentialEquation;

  class MockInstance : public Integrator<ODE>::Instance {
   public:
    MockInstance() : Integrator<ODE>::Instance() {}

    MOCK_METHOD1_T(Solve, Status(Instant const& t_final));

    MOCK_CONST_METHOD1_T(
        WriteToMessage,
        void(not_null<serialization::IntegratorInstance*> message));
  };

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      IntegrationProblem<ODE> const& problem,
      typename Integrator<ODE>::AppendState&& append_state,
      Time const& step) const override {
    return make_not_null_unique<MockInstance>();
  }

  static MockFixedStepSizeIntegrator const& Get() {
    static MockFixedStepSizeIntegrator const integrator;
    return integrator;
  }

 private:
  MockFixedStepSizeIntegrator() : FixedStepSizeIntegrator<ODE>(
      serialization::FixedStepSizeIntegrator::DUMMY) {}
};

}  // namespace internal_ordinary_differential_equations

using internal_ordinary_differential_equations::MockFixedStepSizeIntegrator;

}  // namespace integrators
}  // namespace principia
