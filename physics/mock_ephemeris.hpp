#pragma once

#include <vector>

#include "gmock/gmock.h"
#include "physics/ephemeris.hpp"

namespace principia {
namespace physics {

template<typename Frame>
class MockEphemeris : public Ephemeris<Frame> {
 public:
  MockEphemeris() : Ephemeris<Frame>() {}

  MOCK_CONST_METHOD0_T(bodies, std::vector<MassiveBody const*> const&());
  MOCK_CONST_METHOD1_T(trajectory,
                       not_null<ContinuousTrajectory<Frame> const*>(
                           not_null<MassiveBody const*> body));
  MOCK_CONST_METHOD0_T(empty, bool());
  MOCK_CONST_METHOD0_T(t_min, Instant());
  MOCK_CONST_METHOD0_T(t_max, Instant());
  MOCK_CONST_METHOD0_T(
      planetary_integrator,
      FixedStepSizeIntegrator<
          typename Ephemeris<Frame>::NewtonianMotionEquation> const&());

  MOCK_METHOD1_T(ForgetBefore, void(Instant const& t));
  MOCK_METHOD1_T(Prolong, void(Instant const& t));
  MOCK_METHOD6_T(
      FlowWithAdaptiveStep,
      void(not_null<DiscreteTrajectory<Frame>*> const trajectory,
           typename Ephemeris<Frame>::IntrinsicAcceleration
               intrinsic_acceleration,
           Length const& length_integration_tolerance,
           Speed const& speed_integration_tolerance,
           AdaptiveStepSizeIntegrator<
               typename Ephemeris<Frame>::NewtonianMotionEquation> const&
               integrator,
           Instant const& t));
  MOCK_METHOD4_T(
      FlowWithFixedStep,
      void(std::vector<not_null<DiscreteTrajectory<Frame>*>> const&
               trajectories,
           typename Ephemeris<Frame>::IntrinsicAccelerations const&
               intrinsic_accelerations,
           Time const& step,
           Instant const& t));

  MOCK_CONST_METHOD2_T(
      ComputeGravitationalAcceleration,
      Vector<Acceleration, Frame>(Position<Frame> const& position,
                                  Instant const & t));

  // NOTE(phl): Can't mock the other overloads of
  // ComputeGravitationalAcceleration, it causes an internal error in the
  // compiler.

  MOCK_CONST_METHOD1_T(WriteToMessage,
                       void(not_null<serialization::Ephemeris*> const message));
};

}  // namespace physics
}  // namespace principia
