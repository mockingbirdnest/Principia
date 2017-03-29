
#pragma once

#include <vector>

#include "gmock/gmock.h"
#include "integrators/mock_integrators.hpp"
#include "physics/ephemeris.hpp"

namespace principia {
namespace physics {
namespace internal_ephemeris {

using integrators::MockFixedStepSizeIntegrator;

template<typename Frame>
class MockEphemeris : public Ephemeris<Frame> {
 public:
  using typename Ephemeris<Frame>::AdaptiveStepParameters;
  using typename Ephemeris<Frame>::FixedStepParameters;

  MockEphemeris()
      : Ephemeris<Frame>(
            MockFixedStepSizeIntegrator<
                typename Ephemeris<Frame>::NewtonianMotionEquation>::Get()) {}

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
      bool(not_null<DiscreteTrajectory<Frame>*> trajectory,
           typename Ephemeris<Frame>::IntrinsicAcceleration
               intrinsic_acceleration,
           Instant const& t,
           AdaptiveStepParameters const& parameters,
           std::int64_t max_ephemeris_steps,
           bool last_point_only));
  MOCK_METHOD4_T(
      FlowWithFixedStep,
      void(std::vector<not_null<DiscreteTrajectory<Frame>*>> const&
               trajectories,
           typename Ephemeris<Frame>::IntrinsicAccelerations const&
               intrinsic_accelerations,
           Instant const& t,
           FixedStepParameters const& parameters));

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
                       int(not_null<MassiveBody const*> body));
  MOCK_CONST_METHOD1_T(
      body_for_serialization_index,
      not_null<MassiveBody const*>(int serialization_index));

  MOCK_CONST_METHOD1_T(WriteToMessage,
                       void(not_null<serialization::Ephemeris*> message));
};

}  // namespace internal_ephemeris

using internal_ephemeris::MockEphemeris;

}  // namespace physics

ACTION_P(AppendToDiscreteTrajectories, degrees_of_freedom) {
  for (auto const& trajectory : arg0) {
    trajectory->Append(arg2, degrees_of_freedom);
  }
}

ACTION_P2(AppendToDiscreteTrajectories, time, degrees_of_freedom) {
  for (auto const& trajectory : arg0) {
    trajectory->Append(time, degrees_of_freedom);
  }
}

ACTION_P(AppendToDiscreteTrajectory, degrees_of_freedom) {
  arg0->Append(arg2, degrees_of_freedom);
}

ACTION_P2(AppendToDiscreteTrajectory, time, degrees_of_freedom) {
  arg0->Append(time, degrees_of_freedom);
}

}  // namespace principia
