#pragma once

#include "physics/discrete_trajectory_segment.hpp"

#include <iterator>

#include "astronomy/epoch.hpp"
#include "glog/logging.h"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment {

using astronomy::InfiniteFuture;
using astronomy::InfinitePast;

template<typename Frame>
DiscreteTrajectorySegment<Frame>::DiscreteTrajectorySegment(
    DiscreteTrajectorySegmentIterator<Frame> const self)
    : self_(self) {}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::begin() const {
  return iterator(self_, timeline_.begin());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::end() const {
  // TODO(phl): We probably don't want empty segments.
  if (timeline_.empty()) {
    return iterator(self_, timeline_.begin());
  } else {
    // The decrement/increment ensures that we normalize the end iterator to the
    // next segment or to the end of the trajectory.
    return ++iterator(self_, --timeline_.end());
  }
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reverse_iterator
DiscreteTrajectorySegment<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reverse_iterator
DiscreteTrajectorySegment<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
bool DiscreteTrajectorySegment<Frame>::empty() const {
  return timeline_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectorySegment<Frame>::size() const {
  // NOTE(phl): This assumes that there are no repeated times *within* a
  // segment.  This is enforced by Append.
  return timeline_.size();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::find(Instant const& t) const {
  auto const it = timeline_.find(t);
  if (it == timeline_.end()) {
    return end();
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::lower_bound(Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it == timeline_.end()) {
    return end();
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::upper_bound(Instant const& t) const {
  auto const it = timeline_.upper_bound(t);
  if (it == timeline_.end()) {
    return end();
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
Instant DiscreteTrajectorySegment<Frame>::t_min() const {
  return empty() ? InfiniteFuture : timeline_.cbegin()->first;
}

template<typename Frame>
Instant DiscreteTrajectorySegment<Frame>::t_max() const {
  return empty() ? InfinitePast : timeline_.crbegin()->first;
}

template<typename Frame>
Position<Frame> DiscreteTrajectorySegment<Frame>::EvaluatePosition(
    Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it->first == t) {
    return it->second.position();
  }
  CHECK_LT(t_min(), t);
  CHECK_GT(t_max(), t);
  return GetInterpolation(it).Evaluate(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectorySegment<Frame>::EvaluateVelocity(
    Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it->first == t) {
    return it->second.velocity();
  }
  CHECK_LT(t_min(), t);
  CHECK_GT(t_max(), t);
  return GetInterpolation(it).EvaluateDerivative(t);
}

template<typename Frame>
DegreesOfFreedom<Frame>
DiscreteTrajectorySegment<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it->first == t) {
    return it->second;
  }
  CHECK_LT(t_min(), t);
  CHECK_GT(t_max(), t);
  auto const interpolation = GetInterpolation(it);
  return {interpolation.Evaluate(t), interpolation.EvaluateDerivative(t)};
}

template<typename Frame>
absl::Status DiscreteTrajectorySegment<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  if (!timeline_.empty() && timeline_.cbegin()->first == t) {
    LOG(WARNING) << "Append at existing time " << t << ", time range = ["
                 << timeline_.cbegin()->first << ", "
                 << timeline_.crbegin()->first << "]";
    return absl::OkStatus();
  }
  auto it = timeline_.emplace_hint(timeline_.cend(),
                                   t,
                                   degrees_of_freedom);
  CHECK(++it == timeline_.end())
      << "Append out of order at " << t << ", last time is "
      << timeline_.crbegin()->first;

  // TODO(phl): Downsampling.
  return absl::OkStatus();
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(Instant const& t) {
  ForgetAfter(timeline_.lower_bound(t));
  // TODO(phl): Downsampling.
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(
    typename Timeline::const_iterator const begin) {
  timeline_.erase(begin, timeline_.end());
  // TODO(phl): Downsampling.
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(Instant const& t) {
  ForgetBefore(timeline_.lower_bound(t));
  // TODO(phl): Downsampling.
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(
    typename Timeline::const_iterator const end) {
  timeline_.erase(timeline_.begin(), end);
  // TODO(phl): Downsampling.
}

template<typename Frame>
Hermite3<Instant, Position<Frame>>
DiscreteTrajectorySegment<Frame>::GetInterpolation(
    typename Timeline::const_iterator const& upper) const {
  CHECK(upper != timeline_.cbegin());
  auto const lower = std::prev(upper);
  auto const& [lower_time, lower_degrees_of_freedom] = *lower;
  auto const& [upper_time, upper_degrees_of_freedom] = *upper;
  return Hermite3<Instant, Position<Frame>>{
      {lower_time, upper_time},
      {lower_degrees_of_freedom.position(),
       upper_degrees_of_freedom.position()},
      {lower_degrees_of_freedom.velocity(),
       upper_degrees_of_freedom.velocity()}};
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_begin() const {
  return timeline_.cbegin();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_end() const {
  return timeline_.cend();
}

}  // namespace internal_discrete_trajectory_segment
}  // namespace physics
}  // namespace principia
