#pragma once

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "integrators/integrators.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace _mock_integrators {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_named_quantities;
using namespace principia::integrators::_integrators;
using namespace principia::quantities::_quantities;

template<typename DifferentialEquation>
class MockFixedStepSizeIntegrator
    : public FixedStepSizeIntegrator<DifferentialEquation> {
 public:
  using ODE = DifferentialEquation;
  using typename Integrator<ODE>::AppendState;

  class MockInstance : public Integrator<ODE>::Instance {
   public:
    MockInstance() : Integrator<ODE>::Instance() {}

    MOCK_METHOD(absl::Status, Solve, (Instant const& t_final), (override));
    MOCK_METHOD(not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>,
                Clone,
                (),
                (const, override));
    MOCK_METHOD(FixedStepSizeIntegrator<ODE> const&,
                integrator,
                (),
                (const, override));
  };

  MOCK_METHOD(not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>,
              NewInstance,
              (InitialValueProblem<ODE> const& problem,
               typename Integrator<ODE>::AppendState const& append_state,
               Time const& step),
              (const, override));
  MOCK_METHOD(void,
              WriteToMessage,
              (not_null<serialization::FixedStepSizeIntegrator*> message),
              (const, override));

  static MockFixedStepSizeIntegrator const& Get() {
    static MockFixedStepSizeIntegrator const integrator;
    return integrator;
  }

 private:
  MockFixedStepSizeIntegrator() : FixedStepSizeIntegrator<ODE>() {}
};

}  // namespace internal

using internal::MockFixedStepSizeIntegrator;

}  // namespace _mock_integrators
}  // namespace integrators
}  // namespace principia

namespace principia::integrators {
using namespace principia::integrators::_mock_integrators;
}  // namespace principia::integrators
