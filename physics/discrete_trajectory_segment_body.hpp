#pragma once

#include "glog/logging.h"
#include "physics/discrete_trajectory_segment.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment {

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_begin() const {
  LOG(FATAL) << "NYI";
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_end() const {
  LOG(FATAL) << "NYI";
}

template<typename Frame>
std::int64_t DiscreteTrajectorySegment<Frame>::size() const {
  LOG(FATAL) << "NYI";
}

}  // namespace internal_discrete_trajectory_segment
}  // namespace physics
}  // namespace principia
