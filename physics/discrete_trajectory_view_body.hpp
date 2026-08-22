#pragma once

#include "physics/discrete_trajectory_view.hpp"

#include <iterator>

namespace principia {
namespace physics {
namespace _discrete_trajectory_view {
namespace internal {

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
    : TrajectoryView<Frame>(trajectory,
                            TrajectoryViewTMin(*trajectory, begin, end),
                            TrajectoryViewTMax(*trajectory, begin, end)),
      begin_(begin),
      end_(end) {}

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
  DCHECK(begin_ != end_);
  return *begin_;
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::reference
DiscreteTrajectoryView<Frame>::back() const {
  DCHECK(begin_ != end_);
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
std::int64_t DiscreteTrajectoryView<Frame>::size() const {
  return std::distance(begin_, end_);
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::find(Instant const& t) const {
  if (this->t_min() <= t && t <= this->t_max()) {
    return trajectory().find(t);
  } else {
    return end_;
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::lower_bound(Instant const& t) const {
  // The order of the tests is important to return `end_` for empty views.
  if (this->t_max() < t) {
    return end_;
  } else if (t <= this->t_min()) {
    return begin_;
  } else {
    return trajectory().lower_bound(t);
  }
}

template<typename Frame>
typename DiscreteTrajectoryView<Frame>::iterator
DiscreteTrajectoryView<Frame>::upper_bound(Instant const& t) const {
  // The order of the tests is important to return `end_` for empty views.
  if (this->t_max() <= t) {
    return end_;
  } else if (t < this->t_min()) {
    return begin_;
  } else {
    return trajectory().upper_bound(t);
  }
}

template<typename Frame>
DiscreteTrajectory<Frame> const& DiscreteTrajectoryView<Frame>::trajectory()
    const {
  return dynamic_cast<DiscreteTrajectory<Frame> const&>(
      TrajectoryView<Frame>::trajectory());
}

template<typename Frame>
Instant DiscreteTrajectoryView<Frame>::TrajectoryViewTMin(
    DiscreteTrajectory<Frame> const& trajectory,
    const_iterator const begin,
    const_iterator const end) {
  if (begin == end) {
    return InfiniteFuture;
  } else {
    // The only way that `begin` can be past the end or `end` at the beginning
    // is if they are equal.
    CHECK(begin != trajectory.end());
    CHECK(end != trajectory.begin());
    return begin->time;
  }
}

template<typename Frame>
Instant DiscreteTrajectoryView<Frame>::TrajectoryViewTMax(
    DiscreteTrajectory<Frame> const& trajectory,
    const_iterator const begin,
    const_iterator const end) {
  if (begin == end) {
    return InfinitePast;
  } else {
    // The only way that `begin` can be past the end or `end` at the beginning
    // is if they are equal.
    CHECK(begin != trajectory.end());
    CHECK(end != trajectory.begin());
    return std::prev(end)->time;
  }
}

}  // namespace internal
}  // namespace _discrete_trajectory_view
}  // namespace physics
}  // namespace principia
