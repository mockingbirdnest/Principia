#pragma once

#include "physics/discrete_trajectory2.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory2 {

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::begin() const {
  return (*segments_->begin())->begin();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::end() const {
  return (*segments_->rbegin())->end();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::reverse_iterator
DiscreteTrajectory2<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::reverse_iterator
DiscreteTrajectory2<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::find(Instant const& t) const {
  return iterator();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::lower_bound(Instant const& t) const {
  return iterator();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::upper_bound(Instant const& t) const {
  return iterator();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentRange
DiscreteTrajectory2<Frame>::segments() const {
  return SegmentRange();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::ReverseSegmentRange
DiscreteTrajectory2<Frame>::rsegments() const {
  return ReverseSegmentRange();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::NewSegment() {
  return SegmentIterator();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::DiscreteTrajectory2
DiscreteTrajectory2<Frame>::DetachSegments(iterator begin) {
  return DiscreteTrajectory2();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::AttachSegments(
    DiscreteTrajectory2&& trajectory) {
  return SegmentIterator();
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::DeleteSegments(iterator begin) {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetAfter(Instant const& t) {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetAfter(iterator begin) {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(Instant const& t) {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(iterator end) {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::Append(Instant const& t, DegreesOfFreedom<Frame> const& degrees_of_freedom) {}

template<typename Frame>
Instant DiscreteTrajectory2<Frame>::t_min() const {
  return Instant();
}

template<typename Frame>
Instant DiscreteTrajectory2<Frame>::t_max() const {
  return Instant();
}

template<typename Frame>
Position<Frame> DiscreteTrajectory2<Frame>::EvaluatePosition(Instant const& time) const {
  return Position<Frame>();
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory2<Frame>::EvaluateVelocity(Instant const& time) const {
  return Velocity<Frame>();
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory2<Frame>::EvaluateDegreesOfFreedom(Instant const& time) const {
  return DegreesOfFreedom<Frame>();
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::WriteToMessage(not_null<serialization::DiscreteTrajectory*> message, std::vector<SegmentIterator> const& tracked, std::vector<iterator> const& exact) const {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::WriteToMessage(not_null<serialization::DiscreteTrajectory*> message, iterator begin, iterator end, std::vector<SegmentIterator> const& tracked, std::vector<iterator> const& exact) const {}

template<typename Frame>
template<typename F, >
not_null<std::unique_ptr<DiscreteTrajectory2>> DiscreteTrajectory2<Frame>::ReadFromMessage(serialization::DiscreteTrajectory const& message, std::vector<DiscreteTrajectory2**> const& tracked) {
  return not_null<std::unique_ptr<DiscreteTrajectory2>>();
}

}  // namespace internal_discrete_trajectory2
}  // namespace physics
}  // namespace principia
