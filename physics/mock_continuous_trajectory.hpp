
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

  MOCK_CONST_METHOD1_T(EvaluatePosition, Position<Frame>(Instant const& time));
  MOCK_CONST_METHOD1_T(EvaluateDegreesOfFreedom,
                       DegreesOfFreedom<Frame>(Instant const& time));
};

}  // namespace internal_continuous_trajectory

using internal_continuous_trajectory::MockContinuousTrajectory;

}  // namespace physics
}  // namespace principia
