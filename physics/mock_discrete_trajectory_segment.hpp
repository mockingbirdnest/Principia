#pragma once

#include <cstdint>

#include "gmock/gmock.h"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment.hpp"

namespace principia {
namespace physics {

template<typename Frame>
class MockDiscreteTrajectorySegment : public DiscreteTrajectorySegment<Frame> {
 public:
  MockDiscreteTrajectorySegment() = default;

  MOCK_METHOD(DiscreteTrajectoryIterator<Frame>, begin, (), (const override));
  MOCK_METHOD(DiscreteTrajectoryIterator<Frame>, end, (), (const override));

  MOCK_METHOD(std::int64_t, size, (), (const override));
};

}  // namespace physics
}  // namespace principia
