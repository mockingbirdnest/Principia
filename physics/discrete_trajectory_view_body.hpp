#pragma once

#include "physics/discrete_trajectory_view.hpp"

#include <iterator>

#include "geometry/instant.hpp"

namespace principia {
namespace physics {
namespace _discrete_trajectory_view {
namespace internal {

using namespace principia::geometry::_instant;

template<typename Frame>
DiscreteTrajectoryView<Frame>::DiscreteTrajectoryView(
    not_null<DiscreteTrajectory<Frame> const*> const trajectory)
    : DiscreteTrajectoryView(trajectory,
                             trajectory->begin(),
                             trajectory->end()) {}

template<typename Frame>
DiscreteTrajectoryView<Frame>::DiscreteTrajectoryView(
    not_null<DiscreteTrajectory<Frame> const*> const trajectory,
    const_iterator const begin,
    const_iterator const end)
    : trajectory_(trajectory), begin_(begin), end_(end) {
  if (begin_ == end_) {
    t_min_ = InfiniteFuture;
    t_max_ = InfinitePast;
  } else {
    // The only way that `begin_` can be past the end or `end_` at the beginning
    // is if they are equal.
    CHECK(begin_ != trajectory_->end());
    CHECK(end_ != trajectory_->begin());
    t_min_ = begin->time;
    t_max_ = std::prev(end_)->time;
    CHECK_LE(t_min_, t_max_) << "Overempty ranges not supported";
  }
}

template<typename Frame>
DiscreteTrajectoryView<Frame>::DiscreteTrajectoryView(
    not_null<DiscreteTrajectory<Frame> const*> const trajectory,
    Instant const& t_min,
    Instant const& t_max)
    : DiscreteTrajectoryView(trajectory,
                             trajectory->lower_bound(t_min),
                             trajectory->upper_bound(t_max)) {}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::front() const {
  DCHECK(!empty());
  return *begin_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::back() const {
  DCHECK(!empty());
  return *std::prev(end_);
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::begin() const {
  return begin_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::end() const {
  return end_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reverse_iterator
DiscreteTrajectoryView<Frame>::rbegin() const {
  return reverse_iterator(end_);
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reverse_iterator
DiscreteTrajectoryView<Frame>::rend() const {
  return reverse_iterator(begin_);
}

template<typename Frame>
bool DiscreteTrajectoryView<Frame>::empty() const {
  return t_max_ < t_min_;
}

template<typename Frame>
std::int64_t DiscreteTrajectoryView<Frame>::size() const {
  return std::distance(begin_, end_);
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::find(Instant const& t) const {
  if (t_min_ <= t && t <= t_max_) {
    return trajectory_->find(t);
  } else {
    return end_;
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::lower_bound(Instant const& t) const {
  // The order of the tests is important to return `end_` for empty views.
  if (t_max_ < t) {
    return end_;
  } else if (t <= t_min_) {
    return begin_;
  } else {
    return trajectory_->lower_bound(t);
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::upper_bound(Instant const& t) const {
  // The order of the tests is important to return `end_` for empty views.
  if (t_max_ <= t) {
    return end_;
  } else if (t < t_min_) {
    return begin_;
  } else {
    return trajectory_->upper_bound(t);
  }
}

template<typename Frame>
Instant DiscreteTrajectoryView<Frame>::t_min() const {
  return t_min_;
}

template<typename Frame>
Instant DiscreteTrajectoryView<Frame>::t_max() const {
  return t_max_;
}

template<typename Frame>
Position<Frame> DiscreteTrajectoryView<Frame>::EvaluatePosition(
    Instant const& t) const {
  CHECK_LE(t_min_, t);
  CHECK_LE(t, t_max_);
  return trajectory_->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectoryView<Frame>::EvaluateVelocity(
    Instant const& t) const {
  CHECK_LE(t_min_, t);
  CHECK_LE(t, t_max_);
  return trajectory_->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectoryView<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  CHECK_LE(t_min_, t);
  CHECK_LE(t, t_max_);
  return trajectory_->EvaluateDegreesOfFreedom(t);
}

}  // namespace internal
}  // namespace _discrete_trajectory_view
}  // namespace physics
}  // namespace principia
