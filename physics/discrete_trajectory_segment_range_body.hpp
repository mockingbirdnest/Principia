#pragma once

#include "physics/discrete_trajectory_segment_range.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment_range {

template<typename Iterator>
DiscreteTrajectorySegmentRange<Iterator>::DiscreteTrajectorySegmentRange(
    Iterator const begin,
    Iterator const end)
    : begin_(begin), end_(end) {}

template<typename Iterator>
Iterator DiscreteTrajectorySegmentRange<Iterator>::begin() const {
  return begin_;
}

template<typename Iterator>
Iterator DiscreteTrajectorySegmentRange<Iterator>::end() const {
  return end_;
}

template<typename Iterator>
inline bool DiscreteTrajectorySegmentRange<Iterator>::empty() const {
  return begin_ == end_;
}

template<typename Iterator>
inline std::int64_t DiscreteTrajectorySegmentRange<Iterator>::size() const {
  return std::distance(begin_, end_);
}

}  // namespace internal_discrete_trajectory_segment_range
}  // namespace physics
}  // namespace principia
