
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
  using typename Ephemeris<Frame>::IntrinsicAcceleration;
  using typename Ephemeris<Frame>::IntrinsicAccelerations;
  using typename Ephemeris<Frame>::NewtonianMotionEquation;

  MockEphemeris()
      : Ephemeris<Frame>(
            MockFixedStepSizeIntegrator<NewtonianMotionEquation>::Get()) {}

  MOCK_METHOD(std::vector<not_null<MassiveBody const*>> const&,
              bodies,
              (),
              (const, override));
  MOCK_METHOD(not_null<ContinuousTrajectory<Frame> const*>,
              trajectory,
              (not_null<MassiveBody const*> body),
              (const, override));
  MOCK_METHOD(bool, empty, (), (const, override));
  MOCK_METHOD(Instant, t_min, (), (const, override));
  MOCK_METHOD(Instant, t_max, (), (const, override));
  MOCK_METHOD(FixedStepSizeIntegrator<NewtonianMotionEquation> const&,
              planetary_integrator,
              (),
              (const, override));

  MOCK_METHOD(absl::Status, Prolong, (Instant const& t), (override));
  MOCK_METHOD(
      not_null<std::unique_ptr<
          typename Integrator<NewtonianMotionEquation>::Instance>>,
      NewInstance,
      (std::vector<not_null<DiscreteTraject0ry<Frame>*>> const& trajectories,
       IntrinsicAccelerations const& intrinsic_accelerations,
       FixedStepParameters const& parameters),
      (override));
  MOCK_METHOD(absl::Status,
              FlowWithAdaptiveStep,
              (not_null<DiscreteTraject0ry<Frame>*> trajectory,
               IntrinsicAcceleration intrinsic_acceleration,
               Instant const& t,
               AdaptiveStepParameters const& parameters,
               std::int64_t max_ephemeris_steps),
              (override));
  MOCK_METHOD(
      absl::Status,
      FlowWithFixedStep,
      (Instant const& t,
       typename Integrator<NewtonianMotionEquation>::Instance& instance),
      (override));

  MOCK_METHOD((Vector<Acceleration, Frame>),
              ComputeGravitationalAccelerationOnMasslessBody,
              (Position<Frame> const& position, Instant const& t),
              (const, override));

  MOCK_METHOD((Vector<Acceleration, Frame>),
              ComputeGravitationalAccelerationOnMasslessBody,
              (not_null<DiscreteTraject0ry<Frame>*> trajectory,
               Instant const& t),
              (const, override));

  MOCK_METHOD((Vector<Acceleration, Frame>),
              ComputeGravitationalAccelerationOnMassiveBody,
              (not_null<MassiveBody const*> body,
               Instant const& t),
              (const, override));

  MOCK_METHOD(int,
              serialization_index_for_body,
              (not_null<MassiveBody const*> body),
              (const, override));
  MOCK_METHOD(not_null<MassiveBody const*>,
              body_for_serialization_index,
              (int serialization_index),
              (const, override));

  MOCK_METHOD(void,
              WriteToMessage,
              (not_null<serialization::Ephemeris*> message),
              (const, override));

  MOCK_METHOD(Instant, t_min_locked, (), (const, override));
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

ACTION_P3(AppendToDiscreteTrajectory, trajectory, time, degrees_of_freedom) {
  // The extra level of indirection is useful for tests that get a pointer to a
  // trajectory and squirrel it away using |SaveArg<N>|.
  (*trajectory)->Append(time, degrees_of_freedom);
}

ACTION_P(AppendPointsToDiscreteTrajectory, trajectory) {
  for (auto const& [time, degrees_of_freedom] : *trajectory) {
    arg0->Append(time, degrees_of_freedom);
  }
}

}  // namespace principia
