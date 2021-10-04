#pragma once

#include "physics/discrete_trajectory_segment.hpp"

#include "glog/logging.h"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment {

template<typename Frame>
DiscreteTrajectoryIterator<Frame> DiscreteTrajectorySegment<Frame>::begin()
    const {
  LOG(FATAL) << "NYI";
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame> DiscreteTrajectorySegment<Frame>::end()
    const {
  LOG(FATAL) << "NYI";
}

template<typename Frame>
std::int64_t DiscreteTrajectorySegment<Frame>::size() const {
  LOG(FATAL) << "NYI";
}

}  // namespace internal_discrete_trajectory_segment
}  // namespace physics
}  // namespace principia
