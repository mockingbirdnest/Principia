
#pragma once

#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "integrators/integrators.hpp"

namespace principia {
namespace integrators {
namespace internal_integrators {

using base::make_not_null_unique;

template<typename DifferentialEquation>
class MockFixedStepSizeIntegrator
    : public FixedStepSizeIntegrator<DifferentialEquation> {
 public:
  using ODE = DifferentialEquation;
  using typename Integrator<ODE>::AppendState;

  class MockInstance : public Integrator<ODE>::Instance {
   public:
    MockInstance() : Integrator<ODE>::Instance() {}

    MOCK_METHOD1_T(Solve, Status(Instant const& t_final));
    MOCK_CONST_METHOD0_T(
        Clone,
        not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>());
    MOCK_CONST_METHOD0_T(integrator, FixedStepSizeIntegrator<ODE> const&());
  };

  MOCK_CONST_METHOD3_T(
      NewInstance,
      not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>(
          IntegrationProblem<ODE> const& problem,
          typename Integrator<ODE>::AppendState const& append_state,
          Time const& step));
  MOCK_CONST_METHOD1_T(
      WriteToMessage,
      void(not_null<serialization::FixedStepSizeIntegrator*> message));
  MOCK_CONST_METHOD4_T(
      ReadFromMessage,
      not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>(
          serialization::FixedStepSizeIntegratorInstance const& message,
          IntegrationProblem<ODE> const& problem,
          AppendState const& append_state,
          Time const& step));

  static MockFixedStepSizeIntegrator const& Get() {
    static MockFixedStepSizeIntegrator const integrator;
    return integrator;
  }

 private:
  MockFixedStepSizeIntegrator() : FixedStepSizeIntegrator<ODE>() {}
};

}  // namespace internal_integrators

using internal_integrators::MockFixedStepSizeIntegrator;

}  // namespace integrators
}  // namespace principia
