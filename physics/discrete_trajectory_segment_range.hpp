#pragma once

#include <cstdint>

namespace principia {
namespace physics {
namespace _discrete_trajectory_segment_range {
namespace internal {

// A range of segments in a DiscreteTrajectory, iterator upon using |Iterator|.
// Convenient for range-based loops.
// TODO(phl): Move to base or use the Ranges library if it turns out that this
// class doesn't need to know more about trajectories.
template<typename Iterator>
class DiscreteTrajectorySegmentRange {
 public:
  DiscreteTrajectorySegmentRange() = default;
  DiscreteTrajectorySegmentRange(Iterator begin, Iterator end);

  typename Iterator::reference front() const;
  typename Iterator::reference back() const;

  Iterator begin() const;
  Iterator end() const;

  bool empty() const;
  std::int64_t size() const;

 private:
  Iterator begin_;
  Iterator end_;
};

}  // namespace internal

using internal_discrete_trajectory_segment_range::
      DiscreteTrajectorySegmentRange;

}  // namespace _discrete_trajectory_segment_range
}  // namespace physics
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_discrete_trajectory_segment_range;
}  // namespace principia::physics

#include "physics/discrete_trajectory_segment_range_body.hpp"
