#include "physics/discrete_trajectory_segment_iterator.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment_iterator {

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>&
DiscreteTrajectorySegmentIterator<Frame>::operator++() {
  ++iterator_;
  return *this;
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>&
DiscreteTrajectorySegmentIterator<Frame>::operator--() {
  --iterator_;
  return *this;
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectorySegmentIterator<Frame>::operator++(int) {  // NOLINT
  return DiscreteTrajectorySegmentIterator(segments_, iterator_++);
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectorySegmentIterator<Frame>::operator--(int) {  // NOLINT
  return DiscreteTrajectorySegmentIterator(segments_, iterator_--);
}

template<typename Frame>
internal_discrete_trajectory_segment::DiscreteTrajectorySegment<Frame> const&
DiscreteTrajectorySegmentIterator<Frame>::operator*() const {
  return *iterator_;
}

template<typename Frame>
internal_discrete_trajectory_segment::DiscreteTrajectorySegment<Frame> const*
DiscreteTrajectorySegmentIterator<Frame>::operator->() const {
  return &*iterator_;
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::operator==(
    DiscreteTrajectorySegmentIterator const& other) const {
  return segments_ == other.segments_ && iterator_ == other.iterator_;
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::operator!=(
    DiscreteTrajectorySegmentIterator const& other) const {
  return !operator==(other);
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>::DiscreteTrajectorySegmentIterator(
    not_null<Segments const*> const segments,
    typename Segments::const_iterator iterator)
    : segments_(segments),
      iterator_(iterator) {}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::is_begin() const {
  return iterator_ == segments_->begin();
}

template<typename Frame>
bool DiscreteTrajectorySegmentIterator<Frame>::is_end() const {
  return iterator_ == segments_->end();
}

template<typename Frame>
DiscreteTrajectorySegmentRange<DiscreteTrajectorySegmentIterator<Frame>>
DiscreteTrajectorySegmentIterator<Frame>::segments() const {
  return {DiscreteTrajectorySegmentIterator(segments_, segments_->begin()),
          DiscreteTrajectorySegmentIterator(segments_, segments_->end())};
}

}  // namespace internal_discrete_trajectory_segment_iterator
}  // namespace physics
}  // namespace principia
