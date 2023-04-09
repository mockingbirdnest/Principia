#pragma once

#include "physics/discrete_trajectory_segment_iterator.hpp"

namespace principia {
namespace physics {
namespace _discrete_trajectory_segment_iterator {
namespace internal {

// Note the use of DCHECK, not DCHECK_NOTNULL, below, because the latter does
// not go away when compiled in non-debug mode (don't ask).

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>&
DiscreteTrajectorySegmentIterator<Frame>::operator++() {
  DCHECK(segments_ != nullptr);
  ++iterator_;
  return *this;
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>&
DiscreteTrajectorySegmentIterator<Frame>::operator--() {
  DCHECK(segments_ != nullptr);
  --iterator_;
  return *this;
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectorySegmentIterator<Frame>::operator++(int) {  // NOLINT
  DCHECK(segments_ != nullptr);
  return DiscreteTrajectorySegmentIterator(segments_, iterator_++);
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectorySegmentIterator<Frame>::operator--(int) {  // NOLINT
  DCHECK(segments_ != nullptr);
  return DiscreteTrajectorySegmentIterator(segments_, iterator_--);
}

template<typename Frame>
typename DiscreteTrajectorySegmentIterator<Frame>::reference
DiscreteTrajectorySegmentIterator<Frame>::operator*() const {
  DCHECK(segments_ != nullptr);
  return *iterator_;
}

template<typename Frame>
typename DiscreteTrajectorySegmentIterator<Frame>::pointer
DiscreteTrajectorySegmentIterator<Frame>::operator->() const {
  DCHECK(segments_ != nullptr);
  return &*iterator_;
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::operator==(
    DiscreteTrajectorySegmentIterator const& other) const {
  DCHECK(segments_ != nullptr);
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
  DCHECK(segments_ != nullptr);
  return iterator_ == segments_->begin();
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::is_end() const {
  DCHECK(segments_ != nullptr);
  return iterator_ == segments_->end();
}

template<typename Frame>
DiscreteTrajectorySegmentRange<DiscreteTrajectorySegmentIterator<Frame>>
DiscreteTrajectorySegmentIterator<Frame>::segments() const {
  DCHECK(segments_ != nullptr);
  return {DiscreteTrajectorySegmentIterator(segments_, segments_->begin()),
          DiscreteTrajectorySegmentIterator(segments_, segments_->end())};
}

template<typename Frame>
typename DiscreteTrajectorySegmentIterator<Frame>::Segments::iterator
DiscreteTrajectorySegmentIterator<Frame>::iterator() const {
  return iterator_;
}

}  // namespace internal
}  // namespace _discrete_trajectory_segment_iterator
}  // namespace physics
}  // namespace principia
