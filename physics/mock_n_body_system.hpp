#pragma once

#include "physics/ephemeris.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace physics {

template<typename Frame>
class MockNBodySystem : public Ephemeris<Frame> {
 public:
  MockNBodySystem() : Ephemeris() {}

  MOCK_CONST_METHOD0_T(bodies, std::vector<MassiveBody const*> const&());
  MOCK_CONST_METHOD1_T(
      trajectory,
      not_null<ContinuousTrajectory<Frame> const*>(
          not_null<MassiveBody const*> body));
  MOCK_CONST_METHOD0_T(empty, bool());
  MOCK_CONST_METHOD0_T(t_min, Instant());
  MOCK_CONST_METHOD0_T(t_max, Instant());
  MOCK_CONST_METHOD0_T(
      planetary_integrator,
      FixedStepSizeIntegrator<NewtonianMotionEquation> const&());

  MOCK_METHOD1_T(ForgetBefore, void(Instant const& t));
  MOCK_METHOD1_T(Prolong, void(Instant const& t));
  MOCK_METHOD5_T(FlowWithAdaptiveStep,
                 void(not_null<Trajectory<Frame>*> const trajectory,
                      Length const& length_integration_tolerance,
                      Speed const& speed_integration_tolerance,
                      AdaptiveStepSizeIntegrator<
                          NewtonianMotionEquation> const& integrator,
                      Instant const& t));
  MOCK_METHOD3_T(FlowWithFixedStep,
                 void(Trajectories const& trajectories,
                      Time const& step,
                      Instant const& t));

  MOCK_CONST_METHOD1_T(WriteToMessage,
                       void(not_null<serialization::Ephemeris*> const message));
};

}  // namespace physics
}  // namespace principia
