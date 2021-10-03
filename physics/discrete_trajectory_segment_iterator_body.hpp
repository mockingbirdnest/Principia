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
DiscreteTrajectorySegmentIterator<Frame>::operator++(int) {
  return DiscreteTrajectorySegmentIterator(iterator_++);
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectorySegmentIterator<Frame>::operator--(int) {
  return DiscreteTrajectorySegmentIterator(iterator_--);
}

template<typename Frame>
DiscreteTrajectorySegment<Frame> const&
DiscreteTrajectorySegmentIterator<Frame>::operator*() const {
  return **iterator_;
}

template<typename Frame>
DiscreteTrajectorySegment<Frame> const*
DiscreteTrajectorySegmentIterator<Frame>::operator->() const {
  return iterator_->get();
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>::DiscreteTrajectorySegmentIterator(
    Segments::const_iterator iterator)
    : iterator_(iterator) {}

}  // namespace internal_discrete_trajectory_segment_iterator
}  // namespace physics
}  // namespace principia
