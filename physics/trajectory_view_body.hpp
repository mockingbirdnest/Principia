#pragma once

#include "physics/trajectory_view.hpp"

namespace principia {
namespace physics {
namespace _trajectory_view {
namespace internal {

template<typename Frame>
TrajectoryView<Frame>::TrajectoryView(
    not_null<Trajectory<Frame> const*> const trajectory)
    : TrajectoryView(trajectory, trajectory->t_min(), trajectory->t_max()) {}

template<typename Frame>
TrajectoryView<Frame>::TrajectoryView(
    not_null<Trajectory<Frame> const*> const trajectory,
    Instant const& t_min,
    Instant const& t_max)
    : trajectory_(trajectory), t_min_(t_min), t_max_(t_max) {
  if (t_max_ < t_min_) {
    t_min_ = InfiniteFuture;
    t_max_ = InfinitePast;
  } else {
    CHECK_LE(trajectory_->t_min(), t_min_);
    CHECK_LE(t_max_, trajectory_->t_max());
  }
}

template<typename Frame>
bool TrajectoryView<Frame>::empty() const {
  return t_max_ < t_min_;
}

template<typename Frame>
Instant TrajectoryView<Frame>::t_min() const {
  return t_min_;
}

template<typename Frame>
Instant TrajectoryView<Frame>::t_max() const {
  return t_max_;
}

template<typename Frame>
Position<Frame> TrajectoryView<Frame>::EvaluatePosition(
    Instant const& t) const {
  CHECK_LE(t_min_, t);
  CHECK_LE(t, t_max_);
  return trajectory_->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> TrajectoryView<Frame>::EvaluateVelocity(
    Instant const& t) const {
  CHECK_LE(t_min_, t);
  CHECK_LE(t, t_max_);
  return trajectory_->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> TrajectoryView<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  CHECK_LE(t_min_, t);
  CHECK_LE(t, t_max_);
  return trajectory_->EvaluateDegreesOfFreedom(t);
}

template<typename Frame>
Trajectory<Frame> const& TrajectoryView<Frame>::trajectory() const {
  return *trajectory_;
}

}  // namespace internal
}  // namespace _trajectory_view
}  // namespace physics
}  // namespace principia
