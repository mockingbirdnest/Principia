#include "physics/discrete_trajectory_iterator.hpp"

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_iterator {

using geometry::Instant;

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator++() {
  NormalizeAtSegmentBegin(segment_, point_, previous_time_);
  for (;;) {
    if (IsAtSegmentRBegin()) {
      ++segment_;
      point_ = AtSegmentBegin{};
      break;
    } else {
      auto& point = iterator(point_);
      ++point;
      if (point->first != previous_time_) {
        break;
      }
    }
  }
  return *this;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator--() {
  NormalizeAtSegmentRBegin(segment_, point_, previous_time_);
  for (;;) {
    if (IsAtSegmentBegin()) {
      --segment_;
      point_ = AtSegmentRBegin{};
      break;
    } else {
      auto& point = iterator(point_);
      --point;
      if (point->first != previous_time_) {
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
  auto point = point_;
  Instant time;
  NormalizeAtSegmentTips(segment_, point, time);
  if (time == previous_time_) {
    DiscreteTrajectoryIterator it(segment_, point);
    --it;
    point = it.point_;
    NormalizeAtSegmentTips(it.segment_, point, time);
  }
  return *iterator(point);
}

template<typename Frame>
typename internal_discrete_trajectory_types::Timeline<Frame>::value_type const*
DiscreteTrajectoryIterator<Frame>::operator->() const {
  auto point = point_;
  Instant time;
  NormalizeAtSegmentTips(segment_, point, time);
  if (time == previous_time_) {
    DiscreteTrajectoryIterator it(segment_, point);
    ++it;
    point = it.point_;
    NormalizeAtSegmentTips(it.segment_, point, time);
  }
  return &*iterator(point);
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator==(
    DiscreteTrajectoryIterator const& other) const {
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::operator!=(
    DiscreteTrajectoryIterator const& other) const {
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>::DiscreteTrajectoryIterator(
    DiscreteTrajectorySegmentIterator<Frame> const segment,
    LazyTimelineConstIterator const point)
    : segment_(segment),
      point_(point) {}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::IsAtSegmentBegin() const {
  return std::holds_alternative<AtSegmentBegin>(point_) ||
         (std::holds_alternative<typename Timeline::const_iterator>(point_) &&
          iterator(point_) == segment_->timeline_begin());
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::IsAtSegmentRBegin() const {
  return std::holds_alternative<AtSegmentRBegin>(point_) ||
         (std::holds_alternative<typename Timeline::const_iterator>(point_) &&
          iterator(point_) == --segment_->timeline_end());
}

template<typename Frame>
void DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentBegin(
    DiscreteTrajectorySegmentIterator<Frame> const& segment,
    LazyTimelineConstIterator& point,
    Instant& time) {
  if (std::holds_alternative<AtSegmentBegin>(point)) {
    point = segment->timeline_begin();
  }
  if (!std::holds_alternative<AtSegmentRBegin>(point)) {
    time = iterator(point)->first;
  }
}

template<typename Frame>
void DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentRBegin(
    DiscreteTrajectorySegmentIterator<Frame> const& segment,
    LazyTimelineConstIterator& point,
    Instant& time) {
  if (std::holds_alternative<AtSegmentRBegin>(point)) {
    point = --segment->timeline_end();
  }
  if (!std::holds_alternative<AtSegmentBegin>(point)) {
    time = iterator(point)->first;
  }
}

template<typename Frame>
void DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentTips(
    DiscreteTrajectorySegmentIterator<Frame> const& segment,
    LazyTimelineConstIterator& point,
    Instant& time) {
  if (std::holds_alternative<AtSegmentBegin>(point)) {
    point = segment->timeline_begin();
  }
  if (std::holds_alternative<AtSegmentRBegin>(point)) {
    point = --segment->timeline_end();
  }
  time = iterator(point)->first;
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::Timeline::const_iterator&
DiscreteTrajectoryIterator<Frame>::iterator(
    LazyTimelineConstIterator& point) {
  DCHECK(std::holds_alternative<typename Timeline::const_iterator>(point));
  return std::get<typename Timeline::const_iterator>(point);
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::Timeline::const_iterator const&
DiscreteTrajectoryIterator<Frame>::iterator(
    LazyTimelineConstIterator const& point) {
  DCHECK(std::holds_alternative<typename Timeline::const_iterator>(point));
  return std::get<typename Timeline::const_iterator>(point);
}

}  // namespace internal_discrete_trajectory_iterator
}  // namespace physics
}  // namespace principia
