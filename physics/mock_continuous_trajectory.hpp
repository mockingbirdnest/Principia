#pragma once

#include <vector>

#include "gmock/gmock.h"
#include "physics/continuous_trajectory.hpp"

namespace principia {
namespace physics {

template<typename Frame>
class MockContinuousTrajectory : public ContinuousTrajectory<Frame> {
 public:
  MockContinuousTrajectory() : ContinuousTrajectory<Frame>() {}

  MOCK_CONST_METHOD2_T(
      EvaluateDegreesOfFreedom,
      DegreesOfFreedom<Frame>(
          Instant const& time,
          typename ContinuousTrajectory<Frame>::Hint* const hint));
};

}  // namespace physics
}  // namespace principia
