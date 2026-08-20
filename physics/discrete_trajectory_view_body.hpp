#pragma once

#include "physics/discrete_trajectory_view.hpp"

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_cat.h"
#include "base/status_utilities.hpp"  // 🧙 For RETURN_IF_ERROR.
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace _discrete_trajectory {
namespace internal {

using namespace principia::quantities::_quantities;

template<typename Frame>
DiscreteTrajectoryView<Frame>::DiscreteTrajectoryView(
    DiscreteTrajectory<Frame> const& trajectory,
    const_iterator const begin,
    const_iterator const end)
    : trajectory_(&trajectory), begin_(begin), end_(end) {
  if (begin_ == trajectory_->end()) {
    t_min = InfiniteFuture;
  } else {
    t_min = begin->time;
  }
  if (end_ == trajectory_->begin()) {
    t_max = InfinitePast;
  } else {
    t_max = std::prev(end_)->time;
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::front() const {
  return *begin_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::back() const {
  return *std::prev(end_);
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::begin() const {
  return begin_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::end() const {
  return end_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reverse_iterator
DiscreteTrajectoryView<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reverse_iterator
DiscreteTrajectoryView<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
bool DiscreteTrajectoryView<Frame>::empty() const {
  return begin_ == end_;
}

template<typename Frame>
std::int64_t DiscreteTrajectoryView<Frame>::size() const {
  return std::distance(begin_, end_);
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::find(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    return end();
  }
  auto const sit = leit->second;
  auto const it = sit->FindOrNullopt(t);
  return it.has_value() ? *it : end();
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::lower_bound(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    // This includes an empty trajectory.
    return begin();
  }
  auto const sit = leit->second;
  auto const it = sit->LowerBoundOrNullopt(t);
  return it.has_value() ? *it : end();
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::upper_bound(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    // This includes an empty trajectory.
    return begin();
  }
  auto const sit = leit->second;
  auto const it = sit->UpperBoundOrNullopt(t);
  return it.has_value() ? *it : end();
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::SegmentRange
DiscreteTrajectoryView<Frame>::segments() const {
  return SegmentRange(SegmentIterator(segments_.get(), segments_->begin()),
                      SegmentIterator(segments_.get(), segments_->end()),
                      segments_->size());
}

template<typename Frame>
std::ranges::reverse_view<typename DiscreteTrajectoryView<Frame>::SegmentRange>
DiscreteTrajectoryView<Frame>::rsegments() const {
  return std::ranges::reverse_view(segments());
}

template<typename Frame>
Instant DiscreteTrajectoryView<Frame>::t_min() const {
  return t_min_;
}

template<typename Frame>
Instant DiscreteTrajectoryView<Frame>::t_max() const {
  return t_max_;
}

template<typename Frame>
Position<Frame> DiscreteTrajectoryView<Frame>::EvaluatePosition(
    Instant const& t) const {
  return trajectory_->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectoryView<Frame>::EvaluateVelocity(
    Instant const& t) const {
  return trajectory_->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectoryView<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  return trajectory_->EvaluateDegreesOfFreedom(t);
}

}  // namespace internal
}  // namespace _discrete_trajectory
}  // namespace physics
}  // namespace principia
