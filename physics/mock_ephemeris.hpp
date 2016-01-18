
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

  MOCK_CONST_METHOD0_T(bodies,
                       std::vector<not_null<MassiveBody const*>> const&());
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
      ComputeGravitationalAccelerationOnMasslessBody,
      Vector<Acceleration, Frame>(Position<Frame> const& position,
                                  Instant const & t));

  // NOTE(phl): This overload introduces ambiguities in the expectations.
  // MOCK_CONST_METHOD2_T(
  //     ComputeGravitationalAccelerationOnMasslessBody,
  //     Vector<Acceleration, Frame>(
  //         not_null<DiscreteTrajectory<Frame>*> /*const*/ trajectory,
  //         Instant const& t));

  // NOTE(phl): The commented-out const below is to work-around a compiler
  // internal error.  Don't ask.
  MOCK_CONST_METHOD2_T(
      ComputeGravitationalAccelerationOnMassiveBody,
      Vector<Acceleration, Frame>(
          not_null<MassiveBody const*> /*const*/ body,
          Instant const& t));

  MOCK_CONST_METHOD1_T(serialization_index_for_body,
                       int(not_null<MassiveBody const*> const body));
  MOCK_CONST_METHOD1_T(
      body_for_serialization_index,
      not_null<MassiveBody const*>(int const serialization_index));

  MOCK_CONST_METHOD1_T(WriteToMessage,
                       void(not_null<serialization::Ephemeris*> const message));
};

}  // namespace physics
}  // namespace principia
