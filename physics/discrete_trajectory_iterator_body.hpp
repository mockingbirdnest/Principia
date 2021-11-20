#pragma once

#include "physics/discrete_trajectory_iterator.hpp"

#include "astronomy/epoch.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_iterator {

using astronomy::InfiniteFuture;
using geometry::Instant;

template<typename Frame>
FORCE_INLINE(inline) DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator++() {
  CHECK(!is_at_end(point_));
  auto& point = iterator(point_);
  Instant const previous_time = point->time;
  do {
    ++point;
    if (point == segment_->timeline_end()) {
      do {
        ++segment_;
      } while (!segment_.is_end() && segment_->timeline_empty());

      if (segment_.is_end()) {
        point_.reset();
        break;
      } else {
        point = segment_->timeline_begin();
      }
    }
  } while (point->time == previous_time);
  return *this;
}

template<typename Frame>
FORCE_INLINE(inline) DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator--() {
  bool const point_is_at_end = is_at_end(point_);
  if (point_is_at_end) {
    // Move the iterator to the end of the last segment.
    segment_ = std::prev(segment_.segments().end());
    point_ = segment_->timeline_end();
    // Now proceed with the decrement.
  }
  auto& point = iterator(point_);
  Instant const previous_time = point_is_at_end ? InfiniteFuture : point->time;
  do {
    if (point == segment_->timeline_begin()) {
      CHECK(!segment_.is_begin());
      --segment_;
      point = segment_->timeline_end();
    }
    --point;
  } while (point->time == previous_time);
  return *this;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>
DiscreteTrajectoryIterator<Frame>::operator++(int) {  // NOLINT
  auto const initial = *this;
  ++*this;
  return initial;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>
DiscreteTrajectoryIterator<Frame>::operator--(int) {  // NOLINT
  auto const initial = *this;
  --*this;
  return initial;
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::reference
DiscreteTrajectoryIterator<Frame>::operator*() const {
  CHECK(!is_at_end(point_));
  return *iterator(point_);
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::pointer
DiscreteTrajectoryIterator<Frame>::operator->() const {
  CHECK(!is_at_end(point_));
  return &*iterator(point_);
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator+=(difference_type const n) {
  if (n < 0) {
    return *this -= (-n);
  } else {
    // This loop attempts to skip entire segments.  To do this, it relies on
    // how operator++ moves through the trajectory.  If this was to change this
    // function might become less efficient, but it would not become incorrect
    // (it would fall back to vanilla increments).
    difference_type m = n;
    while (m > 0) {
      CHECK(!is_at_end(point_));
      auto& point = iterator(point_);
      // We know that operator++ never leaves |point_| at |timeline_begin()|.
      // Therefore, to detect that we are in a new segment, we must check for
      // the second point of the segment.
      if (segment_->timeline_size() >= 2 &&
          segment_->timeline_size() <= m + 2 &&
          point == std::next(segment_->timeline_begin())) {
        point = std::prev(segment_->timeline_end());
        m -= segment_->timeline_size() - 2;
      } else {
        ++*this;
        --m;
      }
    }
    return *this;
  }
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator-=(difference_type const n) {
  if (n < 0) {
    return *this += (-n);
  } else {
    difference_type m = n;
    if (m > 0 && is_at_end(point_)) {
      --*this;
      --m;
    }
    // This loop attempts to skip entire segments.  To do this, it relies on
    // how operator-- moves through the trajectory.  If this was to change this
    // function might become less efficient, but it would not become incorrect
    // (it would fall back to vanilla decrements).
    while (m > 0) {
      auto& point = iterator(point_);
      // We know that operator-- never leaves |point_| at
      // |std::prev(timeline_end())|.  Therefore, to detect that we are in a new
      // segment, we must check for the second-to-last point of the segment.
      if (segment_->timeline_size() >= 2 &&
          segment_->timeline_size() <= m + 2 &&
          point == std::prev(std::prev(segment_->timeline_end()))) {
        point = segment_->timeline_begin();
        m -= segment_->timeline_size() - 2;
      } else {
        --*this;
        --m;
      }
    }
    return *this;
  }
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::reference
DiscreteTrajectoryIterator<Frame>::operator[](difference_type const n) const {
  return *(*this + n);
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame> DiscreteTrajectoryIterator<Frame>::operator-(
    typename DiscreteTrajectoryIterator<Frame>::difference_type const n) const {
  auto mutable_it = *this;
  return mutable_it -= n;
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::difference_type
DiscreteTrajectoryIterator<Frame>::operator-(
    DiscreteTrajectoryIterator<Frame> const right) const {
  auto const left = *this;
  auto it = right;
  Instant const left_time =
      is_at_end(left.point_) ? InfiniteFuture : left->time;

  // This code is similar to operator+=.
  difference_type m = 0;
  while (it != left) {
    CHECK(!is_at_end(it.point_));
    auto& point = iterator(it.point_);
    auto const& segment = it.segment_;
    if (segment->timeline_size() >= 2 &&
        std::prev(segment->timeline_end())->time <= left_time &&
        point == std::next(segment->timeline_begin())) {
      point = std::prev(segment->timeline_end());
      m += segment->timeline_size() - 2;
    } else {
      ++it;
      ++m;
    }
  }
  return m;
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator==(
    DiscreteTrajectoryIterator const other) const {
  if (is_at_end(point_)) {
    return segment_ == other.segment_ && is_at_end(other.point_);
  } else if (is_at_end(other.point_)) {
    return false;
  } else {
    return iterator(point_)->time == iterator(other.point_)->time;
  }
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator!=(
    DiscreteTrajectoryIterator const other) const {
  return !operator==(other);
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator<(
    DiscreteTrajectoryIterator const other) const {
  if (is_at_end(point_)) {
    return false;
  } else if (is_at_end(other.point_)) {
    return true;
  } else {
    return iterator(point_)->time < iterator(other.point_)->time;
  }
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator>(
    DiscreteTrajectoryIterator const other) const {
  if (is_at_end(other.point_)) {
    return false;
  } else if (is_at_end(point_)) {
    return true;
  } else {
    return iterator(point_)->time > iterator(other.point_)->time;
  }
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator<=(
    DiscreteTrajectoryIterator const other) const {
  return !operator>(other);
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator>=(
    DiscreteTrajectoryIterator const other) const {
  return !operator<(other);
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>::DiscreteTrajectoryIterator(
    DiscreteTrajectorySegmentIterator<Frame> const segment,
    OptionalTimelineConstIterator const point)
    : segment_(segment),
      point_(point) {
  bool incremented_segment = false;
  while (!segment_.is_end() && segment_->timeline_empty()) {
    ++segment_;
    incremented_segment = true;
  }
  if (segment_.is_end()) {
    point_.reset();
  } else if (incremented_segment) {
    point_ = segment_->timeline_begin();
  }
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::is_at_end(
    OptionalTimelineConstIterator const point) {
  return !point.has_value();
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::Timeline::const_iterator&
DiscreteTrajectoryIterator<Frame>::iterator(
    OptionalTimelineConstIterator& point) {
  DCHECK(point.has_value());
  return point.value();
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::Timeline::const_iterator const&
DiscreteTrajectoryIterator<Frame>::iterator(
    OptionalTimelineConstIterator const& point) {
  DCHECK(point.has_value());
  return point.value();
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame> operator+(
    DiscreteTrajectoryIterator<Frame> const it,
    typename DiscreteTrajectoryIterator<Frame>::difference_type const n) {
  auto mutable_it = it;
  return mutable_it += n;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame> operator+(
    typename DiscreteTrajectoryIterator<Frame>::difference_type const n,
    DiscreteTrajectoryIterator<Frame> const it) {
  auto mutable_it = it;
  return mutable_it += n;
}

}  // namespace internal_discrete_trajectory_iterator
}  // namespace physics
}  // namespace principia
