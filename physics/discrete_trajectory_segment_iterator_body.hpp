#pragma once

#include "physics/discrete_trajectory_segment_iterator.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment_iterator {

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>&
DiscreteTrajectorySegmentIterator<Frame>::operator++() {
  CHECK_NOTNULL(segments_);
  ++iterator_;
  return *this;
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>&
DiscreteTrajectorySegmentIterator<Frame>::operator--() {
  CHECK_NOTNULL(segments_);
  --iterator_;
  return *this;
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectorySegmentIterator<Frame>::operator++(int) {  // NOLINT
  CHECK_NOTNULL(segments_);
  return DiscreteTrajectorySegmentIterator(segments_, iterator_++);
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectorySegmentIterator<Frame>::operator--(int) {  // NOLINT
  CHECK_NOTNULL(segments_);
  return DiscreteTrajectorySegmentIterator(segments_, iterator_--);
}

template<typename Frame>
typename DiscreteTrajectorySegmentIterator<Frame>::reference
DiscreteTrajectorySegmentIterator<Frame>::operator*() const {
  CHECK_NOTNULL(segments_);
  return *iterator_;
}

template<typename Frame>
typename DiscreteTrajectorySegmentIterator<Frame>::pointer
DiscreteTrajectorySegmentIterator<Frame>::operator->() const {
  CHECK_NOTNULL(segments_);
  return &*iterator_;
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::operator==(
    DiscreteTrajectorySegmentIterator const& other) const {
  CHECK_NOTNULL(segments_);
  return segments_ == other.segments_ && iterator_ == other.iterator_;
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::operator!=(
    DiscreteTrajectorySegmentIterator const& other) const {
  return !operator==(other);
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>::DiscreteTrajectorySegmentIterator(
    not_null<Segments*> const segments,
    typename Segments::iterator iterator)
    : segments_(segments),
      iterator_(iterator) {}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::is_begin() const {
  CHECK_NOTNULL(segments_);
  return iterator_ == segments_->begin();
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::is_end() const {
  CHECK_NOTNULL(segments_);
  return iterator_ == segments_->end();
}

template<typename Frame>
DiscreteTrajectorySegmentRange<DiscreteTrajectorySegmentIterator<Frame>>
DiscreteTrajectorySegmentIterator<Frame>::segments() const {
  CHECK_NOTNULL(segments_);
  return {DiscreteTrajectorySegmentIterator(segments_, segments_->begin()),
          DiscreteTrajectorySegmentIterator(segments_, segments_->end())};
}

template<typename Frame>
typename DiscreteTrajectorySegmentIterator<Frame>::Segments::iterator
DiscreteTrajectorySegmentIterator<Frame>::iterator() const {
  return iterator_;
}

}  // namespace internal_discrete_trajectory_segment_iterator
}  // namespace physics
}  // namespace principia
