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
    : segments_(make_not_null_unique<Segments>(1)) {
  auto const sit = segments_->begin();
  auto const self = SegmentIterator(segments_.get(), sit);
  *sit = DiscreteTrajectorySegment<Frame>(self);
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::front() const {
  auto const sit = segment_by_left_endpoint_.begin()->second;
  return *sit->timeline_.begin();
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::back() const {
  return *rbegin();
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::begin() const {
  if (empty()) {
    return end();
  } else {
    auto const sit = segment_by_left_endpoint_.begin()->second;
    return sit->begin();
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::end() const {
  return iterator::EndOfLastSegment(
      SegmentIterator(segments_.get(), segments_->end()));
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
  return segment_by_left_endpoint_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectoryView<Frame>::size() const {
  if (empty()) {
    return 0;
  }
  std::int64_t size = 1;
  std::int64_t nonempty_segments = 0;
  for (auto const& segment : *segments_) {
    if (!segment.empty()) {
      ++nonempty_segments;
      size += segment.size();
    }
  }
  size -= nonempty_segments;  // The junction points.
  return size;
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
  if (empty()) {
    return InfiniteFuture;
  }
  for (auto sit = segments_->begin();; ++sit) {
    if (!sit->empty()) {
      return sit->t_min();
    }
  }
}

template<typename Frame>
Instant DiscreteTrajectoryView<Frame>::t_max() const {
  if (empty()) {
    return InfinitePast;
  }
  return segments_->back().t_max();
}

template<typename Frame>
Position<Frame> DiscreteTrajectoryView<Frame>::EvaluatePosition(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectoryView<Frame>::EvaluateVelocity(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectoryView<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluateDegreesOfFreedom(t);
}

}  // namespace internal
}  // namespace _discrete_trajectory
}  // namespace physics
}  // namespace principia
