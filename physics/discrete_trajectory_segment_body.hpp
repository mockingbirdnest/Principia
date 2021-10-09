#pragma once

#include "physics/discrete_trajectory_segment.hpp"

#include "glog/logging.h"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment {

template<typename Frame>
DiscreteTrajectorySegment<Frame>::DiscreteTrajectorySegment(
    DiscreteTrajectorySegmentIterator<Frame> const self)
    : self_(self) {}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::begin() const {
  return iterator(self_, timeline_.begin());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::end() const {
  // TODO(phl): or begin of next segment?
  return iterator(self_, timeline_.end());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reverse_iterator
DiscreteTrajectorySegment<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reverse_iterator
DiscreteTrajectorySegment<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::find(Instant const& t) const {
  return iterator(self_, timeline_.find(t));
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::lower_bound(Instant const& t) const {
  return iterator(self_, timeline_.lower_bound(t));
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::upper_bound(Instant const& t) const {
  return iterator(self_, timeline_.upper_bound(t));
}

template<typename Frame>
bool DiscreteTrajectorySegment<Frame>::empty() const {
  return timeline_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectorySegment<Frame>::size() const {
  // NOTE(phl): This assumes that there are no repeated times *within* a
  // segment.  This is enforced by Append.
  return timeline_.size();
}

template<typename Frame>
absl::Status DiscreteTrajectorySegment<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  if (!timeline_.empty() && timeline_.cbegin()->first == t) {
    LOG(WARNING) << "Append at existing time " << t << ", time range = ["
                 << timeline_.cbegin()->first << ", "
                 << timeline_.crbegin()->first << "]";
    return absl::OkStatus();
  }
  auto it = timeline_.emplace_hint(timeline_.cend(),
                                   t,
                                   degrees_of_freedom);
  CHECK(++it == timeline_.end())
      << "Append out of order at " << t << ", last time is "
      << timeline_.crbegin()->first;

  // TODO(phl): Downsampling.
  return absl::OkStatus();
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(Instant const& t) {}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(
    typename Timeline::const_iterator const begin) {}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(Instant const& t) {}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(
    typename Timeline::const_iterator const end) {}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_begin() const {
  return timeline_.cbegin();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_end() const {
  return timeline_.cend();
}

}  // namespace internal_discrete_trajectory_segment
}  // namespace physics
}  // namespace principia
