#pragma once

#include "physics/discrete_trajectory2.hpp"

#include <vector>

#include "astronomy/epoch.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory2 {

using astronomy::InfinitePast;
using base::make_not_null_unique;
using base::uninitialized;

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
  return ReverseSegmentRange(std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->end())),
                             std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->begin())));
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
DiscreteTrajectory2<Frame>::DetachSegments(SegmentIterator const begin) {
  DiscreteTrajectory2 detached(uninitialized);

  // Move the detached segments to the new trajectory.
  detached.segments_->splice(detached.segments_->end(),
                             *segments_,
                             begin.iterator(), segments_->end());

  // Iterate through the detached segments to move the time-to-segment mapping
  // to the detached trajectory.
  auto endpoint_it = segment_by_left_endpoint_.end();
  for (auto sit = detached.segments_->rbegin();
       sit != detached.segments_->rend();
       ++sit) {
    --endpoint_it;
    detached.segment_by_left_endpoint_.insert(
        segment_by_left_endpoint_.extract(endpoint_it));
  }

  // Reset the self pointers of the new segments.
  for (auto sit = detached.segments_->begin();
       sit != detached.segments_->end();
       ++sit) {
    sit->SetSelf(SegmentIterator(detached.segments_.get(), sit));
  }

  return detached;
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::AttachSegments(
    DiscreteTrajectory2&& trajectory) {
  CHECK(!trajectory.empty());
  // NOTE(phl): This check might be too strict, we might want to allow LT as the
  // time comparison, and to adjust the first point of trajectory as needed.
  // We'll see if the clients need that.
  CHECK_EQ(rbegin()->first, trajectory.begin()->first)
      << "Mismatching times when attaching segments";
  CHECK_EQ(rbegin()->second, trajectory.begin()->second)
      << "Mismatching degrees of freedom when attaching segments";

  if (empty()) {
    *this = DiscreteTrajectory2(uninitialized);
  }

  // The |end| iterator keeps pointing at the end after the splice.  Instead,
  // we track the iterator to the last segment.
  auto const last_before_splice = --segments_->end();

  // Move the attached segments to the end of this trajectory.
  segments_->splice(segments_->end(),
                    *trajectory.segments_);

  auto const end_before_splice = std::next(last_before_splice);
  auto const rbegin_before_splice = std::reverse_iterator(end_before_splice);

  auto endpoint_it = trajectory.segment_by_left_endpoint_.end();
  for (auto sit = segments_->rbegin(); sit != rbegin_before_splice; ++sit) {
    --endpoint_it;
    segment_by_left_endpoint_.insert(
        trajectory.segment_by_left_endpoint_.extract(endpoint_it));
  }

  // Reset the self pointers of the new segments.
  for (auto sit = end_before_splice; sit != segments_->end(); ++sit) {
    sit->SetSelf(SegmentIterator(segments_.get(), sit));
  }

  return SegmentIterator(segments_.get(), end_before_splice);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::DeleteSegments(SegmentIterator const begin) {
  segments_->erase(begin.iterator(), segments_->end());
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetAfter(Instant const& t) {
  auto const sit = FindSegment(t);
  sit->ForgetAfter(t);
  // Here |sit| designates a segment starting at or after |t|.  If |t| is
  // exactly at the beginning of the segment, |ForgetAfter| above will leave it
  // empty.  In that case we drop the segment entirely.  Note that this
  // situation doesn't arise for |ForgetBefore| because of the way |FindSegment|
  // works.
  if (sit->empty()) {
    segments_->erase(sit, segments_->end());
  } else {
    segments_->erase(std::next(sit), segments_->end());
  }
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(Instant const& t) {
  auto const sit = FindSegment(t);
  sit->ForgetBefore(t);
  segments_->erase(segments_->begin(), sit);
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
DiscreteTrajectory2<Frame>::DiscreteTrajectory2(uninitialized_t)
    : segments_(make_not_null_unique<Segments>()) {}

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
