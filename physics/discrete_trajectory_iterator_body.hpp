#include "physics/discrete_trajectory_iterator.hpp"

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_iterator {

using geometry::Instant;

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator++() {
  point_ = NormalizeAtSegmentBegin(point_);
  auto& point = iterator(point_);
  for (;;) {
    ++point;
    if (point == iterator(segment_->end().point_)) {
      ++segment_;
      point_ = AtSegmentBegin{};
      break;
    } else if (point->first != previous_time_) {
      previous_time_ = point->first;
      break;
    }
  }
  return *this;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator--() {
  point_ = NormalizeAtSegmentRBegin(point_);
  auto& point = iterator(point_);
  for (;;) {
    if (point == iterator(segment_->begin().point_)) {
      --segment_;
      point_ = AtSegmentRBegin{};
      break;
    } else {
      --point;
      if (point->first != previous_time_) {
        previous_time_ = point->first;
        break;
      }
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
  //TODO(phl):Do both
  auto const point = NormalizeAtSegmentBegin(NormalizeAtSegmentRBegin(point_));
  return *iterator(point);
}

template<typename Frame>
typename internal_discrete_trajectory_types::Timeline<Frame>::value_type const*
DiscreteTrajectoryIterator<Frame>::operator->() const {
  auto const point = NormalizeAtSegmentBegin(NormalizeAtSegmentRBegin(point_));
  return &*iterator(point);
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>::DiscreteTrajectoryIterator(
    DiscreteTrajectorySegmentIterator<Frame> const segment,
    typename Timeline::const_iterator const point)
    : segment_(segment),
      point_(point) {}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::LazyTimelineConstIterator
DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentBegin(
    LazyTimelineConstIterator const point) const {
  DCHECK(!std::holds_alternative<AtSegmentRBegin>(point));///???
  if (std::holds_alternative<AtSegmentBegin>(point)) {
    return iterator(segment_->begin().point_);
  } else {
    return point;
  }
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::LazyTimelineConstIterator
DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentRBegin(
    LazyTimelineConstIterator const point) const {
  DCHECK(!std::holds_alternative<AtSegmentBegin>(point));///???
  if (std::holds_alternative<AtSegmentRBegin>(point)) {
    auto end = iterator(segment_->end().point_);
    return --end;
  } else {
    return point;
  }
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::Timeline::const_iterator&
DiscreteTrajectoryIterator<Frame>::iterator(
    LazyTimelineConstIterator& point) {
  return std::get<typename Timeline::const_iterator>(point);
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::Timeline::const_iterator const&
DiscreteTrajectoryIterator<Frame>::iterator(
    LazyTimelineConstIterator const& point) {
  return std::get<typename Timeline::const_iterator>(point);
}

}  // namespace internal_discrete_trajectory_iterator
}  // namespace physics
}  // namespace principia
