#pragma once

#include "physics/discrete_trajectory_segment.hpp"

#include <algorithm>
#include <iterator>
#include <list>
#include <vector>

#include "astronomy/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "numerics/fit_hermite_spline.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment {

using astronomy::InfiniteFuture;
using astronomy::InfinitePast;
using geometry::Position;
using numerics::FitHermiteSpline;

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

  if (downsampling_parameters_.has_value()) {
    return DownsampleIfNeeded();
  } else {
    return absl::OkStatus();
  }
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(Instant const& t) {
  ForgetAfter(timeline_.lower_bound(t));
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(
    typename Timeline::const_iterator const begin) {
  std::int64_t number_of_points_to_remove =
      std::distance(begin, timeline_.cend());
  number_of_dense_points_ =
      std::max(0LL, number_of_dense_points_ - number_of_points_to_remove);

  timeline_.erase(begin, timeline_.cend());
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(Instant const& t) {
  ForgetBefore(timeline_.lower_bound(t));
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(
    typename Timeline::const_iterator const end) {
  std::int64_t number_of_points_to_remove =
      std::distance(timeline_.cbegin(), end);
  number_of_dense_points_ =
      std::max(0LL, number_of_dense_points_ - number_of_points_to_remove);

  timeline_.erase(timeline_.cbegin(), end);
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetDownsampling(
    DownsamplingParameters const& downsampling_parameters) {
  // TODO(phl): Do we need this precondition?
  CHECK(!downsampling_parameters_.has_value());
  downsampling_parameters_ = downsampling_parameters;
  number_of_dense_points_ = timeline_.empty() ? 0 : 1;
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ClearDownsampling() {
  downsampling_parameters_ = std::nullopt;
}

template<typename Frame>
absl::Status DiscreteTrajectorySegment<Frame>::DownsampleIfNeeded() {
  ++number_of_dense_points_;
  // Points, hence one more than intervals.
  if (number_of_dense_points_ >
      downsampling_parameters_->max_dense_intervals) {
    // Obtain iterators for all the dense points of the segment.
    using ConstIterators = std::vector<typename Timeline::const_iterator>;
    ConstIterators dense_iterators(number_of_dense_points_);
    CHECK_LE(dense_iterators.size(), timeline_.size());
    auto it = timeline_.crbegin();
    for (int i = dense_iterators.size() - 1; i >= 0; --i) {
      dense_iterators[i] = std::prev(it.base());
      ++it;
    }

    absl::StatusOr<std::list<ConstIterators::const_iterator>>
        right_endpoints = FitHermiteSpline<Instant, Position<Frame>>(
            dense_iterators,
            [](auto&& it) -> auto&& { return it->first; },
            [](auto&& it) -> auto&& { return it->second.position(); },
            [](auto&& it) -> auto&& { return it->second.velocity(); },
            downsampling_parameters_->tolerance);
    if (!right_endpoints.ok()) {
      // Note that the actual appending took place; the propagated status only
      // reflects a lack of downsampling.
      return right_endpoints.status();
    }

    if (right_endpoints->empty()) {
      number_of_dense_points_ = timeline_.empty() ? 0 : 1;
      return absl::OkStatus();
    }

    // Obtain the times for the right endpoints.  This is necessary because we
    // cannot use iterators for erasing points, as they would get invalidated
    // after the first erasure.
    std::vector<Instant> right_endpoints_times;
    right_endpoints_times.reserve(right_endpoints.value().size());
    for (auto const& it_in_dense_iterators : right_endpoints.value()) {
      right_endpoints_times.push_back((*it_in_dense_iterators)->first);
    }

    // Poke holes in the timeline at the places given by
    // |right_endpoints_times|.  This requires one lookup per erasure.
    auto left_it = dense_iterators.front();
    for (Instant const& right : right_endpoints_times) {
      ++left_it;
      auto const right_it = timeline_.find(right);
      left_it = timeline_.erase(left_it, right_it);
    }
    number_of_dense_points_ = std::distance(left_it, timeline_.cend());
  }
  return absl::OkStatus();
}

template<typename Frame>
Hermite3<Instant, Position<Frame>>
DiscreteTrajectorySegment<Frame>::GetInterpolation(
    typename Timeline::const_iterator const upper) const {
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
