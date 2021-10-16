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
  auto const it = segments_->begin();
  *it = DiscreteTrajectorySegment<Frame>(
      DiscreteTrajectorySegmentIterator<Frame>(segments_.get(), it));
  segment_by_left_endpoint_.emplace(InfinitePast, it);
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
  return std::int64_t();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::find(Instant const& t) const {
  return FindSegment(t).find(t);
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::lower_bound(Instant const& t) const {
  return FindSegment(t).lower_bound(t);
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::upper_bound(Instant const& t) const {
  return FindSegment(t).upper_bound(t);
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentRange
DiscreteTrajectory2<Frame>::segments() const {
  return SegmentRange(segments_->begin(), segments_->end());
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::ReverseSegmentRange
DiscreteTrajectory2<Frame>::rsegments() const {
  return ReverseSegmentRange(segments_->rbegin(), segments_->rend());
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::NewSegment() {
  auto& last_segment = segments_->back();
  CHECK(!last_segment.empty())
      << "Cannot create a new segment after an empty one";

  auto const& new_segment = segments_->emplace_back();
  auto const new_segment_it = --segments_->end();
  *new_segment_it = DiscreteTrajectorySegment<Frame>(
      DiscreteTrajectorySegmentIterator<Frame>(segments_.get(),
                                               new_segment_it));

  auto const& [last_time, last_degrees_of_freedom] = *last_segment.rbegin();
  new_segment_it->Append(last_time, last_degrees_of_freedom);
  segment_by_left_endpoint_.emplace_hint(
      segment_by_left_endpoint_.end(), last_time, new_segment_it);

  return SegmentIterator(segments_.get(), new_segment_it);
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
void DiscreteTrajectory2<Frame>::ForgetAfter(Instant const& t) {
  return FindSegment(t).ForgetAfter(t);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetAfter(iterator begin) {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(Instant const& t) {
  return FindSegment(t).ForgetBefore(t);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(iterator end) {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  FindSegment(t).Append(t, degrees_of_freedom);
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
  return FindSegment(t).EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory2<Frame>::EvaluateVelocity(
    Instant const& t) const {
  return FindSegment(t).EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory2<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  return FindSegment(t).EvaluateDegreesOfFreedom(t);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> message,
    std::vector<SegmentIterator> const& tracked,
    std::vector<iterator> const& exact) const {}

template<typename Frame>
void DiscreteTrajectory2<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> message,
    iterator begin,
    iterator end,
    std::vector<SegmentIterator> const& tracked,
    std::vector<iterator> const& exact) const {}

template<typename Frame>
template<typename F, typename>
not_null<std::unique_ptr<DiscreteTrajectory2<Frame>>>
DiscreteTrajectory2<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory2**> const& tracked) {
  return not_null<std::unique_ptr<DiscreteTrajectory2>>();
}

template<typename Frame>
DiscreteTrajectorySegment<Frame>& DiscreteTrajectory2<Frame>::FindSegment(
    Instant const& t) {
  auto const it = --segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin() ||
        segment_by_left_endpoint_.size() == 1)
      << "Time " << t << " in the sentinel segment out of "
      << segment_by_left_endpoint_.size() << " segments";
  return *(it->second);
}

template<typename Frame>
DiscreteTrajectorySegment<Frame> const& DiscreteTrajectory2<Frame>::FindSegment(
    Instant const& t) const {
  auto const it = --segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin() ||
        segment_by_left_endpoint_.size() == 1)
      << "Time " << t << " in the sentinel segment out of "
      << segment_by_left_endpoint_.size() << " segments";
  return *(it->second);
}

}  // namespace internal_discrete_trajectory2
}  // namespace physics
}  // namespace principia
