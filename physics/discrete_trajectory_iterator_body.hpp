#include "physics/discrete_trajectory_iterator.hpp"

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_iterator {

using geometry::Instant;

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator++() {
  Instant const previous_time = point_->first;
  for (;;) {
    ++point_;
    if (point_ == segment_->end().point_) {
      ++segment_;
      //TODO(phl): segment is end?
      point_ = segment_->begin().point_;

      // Skip duplicated points at the beginning of the timeline.
      if (point_->first != previous_time) {
        break;
      }
    } else {
      break;
    }
  }
  return *this;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator--() {
  Instant const previous_time = point_->first;
  for (;;) {
    if (point_ == segment_->begin().point_) {
      //TODO(phl): segment is begin?
      --segment_;
      point_ = --segment_->end().point_;

      // Skip duplicated points at the end of the timeline.
      if (point_->first != previous_time) {
        break;
      }
    } else {
      --point_;
      break;
    }
  }
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
