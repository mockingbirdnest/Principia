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
  if (begin_ == trajectory_->end()) {
    t_min_ = InfiniteFuture;
  } else {
    t_min_ = begin->time;
  }
  if (end_ == trajectory_->begin()) {
    t_max_ = InfinitePast;
  } else {
    t_max_ = std::prev(end_)->time;
  }
  if (t_max_ < t_min_) {
    t_min_ = InfiniteFuture;
    t_max_ = InfinitePast;
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
  return *begin_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::back() const {
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
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reverse_iterator
DiscreteTrajectoryView<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
bool DiscreteTrajectoryView<Frame>::empty() const {
  return t_max_ < t_min_;
}

template<typename Frame>
std::int64_t DiscreteTrajectoryView<Frame>::size() const {
  if (empty()) {
    return 0;
  } else {
    return std::distance(begin(), end());
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::find(Instant const& t) const {
  if (t_min_ <= t && t <= t_max_) {
    return trajectory_->find(t);
  } else {
    return end();
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::lower_bound(Instant const& t) const {
  if (t <= t_max_) {
    return trajectory_->lower_bound(t);
  } else {
    return end();
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::upper_bound(Instant const& t) const {
  if (t < t_max_) {
    return trajectory_->upper_bound(t);
  } else {
    return end();
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
