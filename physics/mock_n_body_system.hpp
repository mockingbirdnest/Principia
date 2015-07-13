#pragma once

#include "physics/ephemeris.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace physics {

template<typename InertialFrame>
class MockNBodySystem : public Ephemeris<InertialFrame> {
  class DummyIntegrator
      : public FixedStepSizeIntegrator<NewtonianMotionEquation> {
    DummyIntegrator() = default;

   public:
    void Solve(IntegrationProblem<ODE> const& problem,
               Time const& step) const override {
      LOG(FATAL) << "dummy";
    }

    static DummyIntegrator const& Instance() {
      static DummyIntegrator const instance;
      return instance;
    }
  };
 public:
  MockNBodySystem() : planetary_integrator_(DummyIntegrator::Instance()) {}

  // TODO(egg): MOCK_ALL_THE_THINGS
};

}  // namespace physics
}  // namespace principia
