#pragma once

#include "physics/ephemeris.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace physics {

template<typename InertialFrame>
class MockNBodySystem : public Ephemeris<InertialFrame> {
 public:
  MockNBodySystem() : Ephemeris() {}

  MOCK_CONST_METHOD0(bodies, std::vector<MassiveBody const*> const&());
  MOCK_CONST_METHOD1(
      trajectory,
      not_null<ContinuousTrajectory<Frame> const*>(
          not_null<MassiveBody const*> body));
  MOCK_CONST_METHOD0(empty, bool());
  MOCK_CONST_METHOD0(t_min, Instant());
  MOCK_CONST_METHOD0(t_max, Instant());
  MOCK_CONST_METHOD0(planetary_integrator,
                     FixedStepSizeIntegrator<NewtonianMotionEquation> const&());

  MOCK_METHOD_1(ForgetBefore, void(Instant const& t));
  MOCK_METHOD1(Prolong, void(Instant const& t));
  MOCK_METHOD5(FlowWithAdaptiveStep,
               void(not_null<Trajectory<Frame>*> const trajectory,
                    Length const& length_integration_tolerance,
                    Speed const& speed_integration_tolerance,
                    AdaptiveStepSizeIntegrator<
                        NewtonianMotionEquation> const& integrator,
                    Instant const& t));
  MOCK_METHOD3(FlowWithFixedStep,
               void(Trajectories const& trajectories,
                    Time const& step,
                    Instant const& t));

  MOCK_CONST_METHOD1(WriteToMessage,
                     void(not_null<serialization::Ephemeris*> const message));
};

}  // namespace physics
}  // namespace principia
