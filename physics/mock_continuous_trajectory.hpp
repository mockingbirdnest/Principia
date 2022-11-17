#pragma once

#include <vector>

#include "gmock/gmock.h"
#include "physics/continuous_trajectory.hpp"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

template<typename Frame>
class MockContinuousTrajectory : public ContinuousTrajectory<Frame> {
 public:
  MockContinuousTrajectory() : ContinuousTrajectory<Frame>() {}

  MOCK_METHOD(Position<Frame>,
              EvaluatePosition,
              (Instant const& time),
              (const, override));
  MOCK_METHOD(DegreesOfFreedom<Frame>,
              EvaluateDegreesOfFreedom,
              (Instant const& time),
              (const, override));
};

}  // namespace internal_continuous_trajectory

using internal_continuous_trajectory::MockContinuousTrajectory;

}  // namespace physics
}  // namespace principia
