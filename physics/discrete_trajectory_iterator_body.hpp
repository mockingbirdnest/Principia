#include "physics/discrete_trajectory_iterator.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_iterator {

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator++() {
  ++point_;
  if (point_ == segment_->end().point_) {
    ++segment_;
    point_ = segment_->begin().point_;
  }
  return *this;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator--() {
  if (point_ == segment_->begin().point_) {
    --segment_;
    point_ = segment_->end().point_;
  }
  --point_;
  return *this;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>
DiscreteTrajectoryIterator<Frame>::operator++(int) {
  auto const initial = *this;
  ++*this;
  return initial;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>
DiscreteTrajectoryIterator<Frame>::operator--(int) {
  auto const initial = *this;
  --*this;
  return initial;
}

template<typename Frame>
typename internal_discrete_trajectory_types::Timeline<Frame>::value_type const&
DiscreteTrajectoryIterator<Frame>::operator*() const {
  return *point_;
}

template<typename Frame>
typename internal_discrete_trajectory_types::Timeline<Frame>::value_type const*
DiscreteTrajectoryIterator<Frame>::operator->() const {
  return &*point_;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>::DiscreteTrajectoryIterator(
    DiscreteTrajectorySegmentIterator<Frame> const segment,
    typename Timeline::const_iterator const point)
    : segment_(segment),
      point_(point) {}

}  // namespace internal_discrete_trajectory_iterator
}  // namespace physics
}  // namespace principia
