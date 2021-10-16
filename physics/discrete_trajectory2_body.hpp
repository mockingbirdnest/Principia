#pragma once

#include "astronomy/epoch.hpp"
#include "physics/discrete_trajectory2.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory2 {

using astronomy::InfinitePast;
using base::make_not_null_unique;

template<typename Frame>
DiscreteTrajectory2<Frame>::DiscreteTrajectory2()
    : segments_(make_not_null_unique<Segments>(1)) {
  auto const sit = segments_->begin();
  *sit =
      DiscreteTrajectorySegment<Frame>(SegmentIterator(segments_.get(), sit));
  segment_by_left_endpoint_.emplace(InfinitePast, sit);
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::begin() const {
  return segments_->front().begin();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::end() const {
  return segments_->back().end();
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
bool DiscreteTrajectory2<Frame>::empty() const {
  for (auto const& segment : *segments_) {
    if (!segment.empty()) {
      return false;
    }
  }
  return true;
}

template<typename Frame>
std::int64_t DiscreteTrajectory2<Frame>::size() const {
  std::int64_t size = 1;
  for (auto const& segment : *segments_) {
    size += segment.size();
  }
  size -= segments_->size();  // The junction points.
  return size;
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::find(Instant const& t) const {
  auto const sit = FindSegment(t);
  auto const it = sit->find(t);
  if (it == sit->end()) {
    return end();
  } else {
    return it;
  }
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::lower_bound(Instant const& t) const {
  auto const sit = FindSegment(t);
  auto const it = sit->lower_bound(t);
  if (it == sit->end()) {
    return end();
  } else {
    return it;
  }
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::upper_bound(Instant const& t) const {
  auto const sit = FindSegment(t);
  auto const it = sit->upper_bound(t);
  if (it == sit->end()) {
    return end();
  } else {
    return it;
  }
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentRange
DiscreteTrajectory2<Frame>::segments() const {
  return SegmentRange(SegmentIterator(
                          segments_.get(), segments_->begin()),
                      SegmentIterator(
                          segments_.get(), segments_->end()));
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::ReverseSegmentRange
DiscreteTrajectory2<Frame>::rsegments() const {
  // TODO(phl): Implement.
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::NewSegment() {
  auto& last_segment = segments_->back();
  CHECK(!last_segment.empty())
      << "Cannot create a new segment after an empty one";

  auto const& new_segment = segments_->emplace_back();
  auto const new_segment_sit = --segments_->end();
  *new_segment_sit = DiscreteTrajectorySegment<Frame>(
      SegmentIterator(segments_.get(), new_segment_sit));

  auto const& [last_time, last_degrees_of_freedom] = *last_segment.rbegin();
  new_segment_sit->Append(last_time, last_degrees_of_freedom);
  segment_by_left_endpoint_.emplace_hint(
      segment_by_left_endpoint_.end(), last_time, new_segment_sit);

  return SegmentIterator(segments_.get(), new_segment_sit);
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::DiscreteTrajectory2
DiscreteTrajectory2<Frame>::DetachSegments(iterator begin) {
  // TODO(phl): Implement.
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::AttachSegments(
    DiscreteTrajectory2&& trajectory) {
  // TODO(phl): Implement.
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::DeleteSegments(iterator begin) {
  // TODO(phl): Implement.
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetAfter(Instant const& t) {
  // TODO(phl): Drop segments as needed.
  return FindSegment(t)->ForgetAfter(t);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetAfter(iterator begin) {
  // TODO(phl): Implement.
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(Instant const& t) {
  // TODO(phl): Drop segments as needed.
  return FindSegment(t)->ForgetBefore(t);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(iterator end) {
  // TODO(phl): Implement.
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  auto sit = FindSegment(t);
  // If this is the first point appended to this trajectory, insert a proper
  // left endpoint and remove the sentinel.
  if (empty()) {
    segment_by_left_endpoint_.emplace_hint(
        segment_by_left_endpoint_.end(), t, sit);
    segment_by_left_endpoint_.erase(segment_by_left_endpoint_.begin());
  }
  sit->Append(t, degrees_of_freedom);
}

template<typename Frame>
Instant DiscreteTrajectory2<Frame>::t_min() const {
  return segments_->front().t_min();
}

template<typename Frame>
Instant DiscreteTrajectory2<Frame>::t_max() const {
  return segments_->back().t_max();
}

template<typename Frame>
Position<Frame> DiscreteTrajectory2<Frame>::EvaluatePosition(
    Instant const& t) const {
  return FindSegment(t)->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory2<Frame>::EvaluateVelocity(
    Instant const& t) const {
  return FindSegment(t)->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory2<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  return FindSegment(t)->EvaluateDegreesOfFreedom(t);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> message,
    std::vector<SegmentIterator> const& tracked,
    std::vector<iterator> const& exact) const {
  // TODO(phl): Implement.
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> message,
    iterator begin,
    iterator end,
    std::vector<SegmentIterator> const& tracked,
    std::vector<iterator> const& exact) const {
  // TODO(phl): Implement.
}

template<typename Frame>
template<typename F, typename>
not_null<std::unique_ptr<DiscreteTrajectory2<Frame>>>
DiscreteTrajectory2<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory2**> const& tracked) {
  // TODO(phl): Implement.
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::Segments::iterator
DiscreteTrajectory2<Frame>::FindSegment(
    Instant const& t) {
  auto it = segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin()) << "No segment covering " << t;
  return (--it)->second;
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::Segments::const_iterator
DiscreteTrajectory2<Frame>::FindSegment(
    Instant const& t) const {
  auto it = segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin()) << "No segment covering " << t;
  return (--it)->second;
}

}  // namespace internal_discrete_trajectory2
}  // namespace physics
}  // namespace principia
