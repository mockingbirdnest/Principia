#include "physics/discrete_trajectory_iterator.hpp"

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_iterator {

using geometry::Instant;

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator++() {
  NormalizeAtSegmentBegin(point_, previous_time_);
  auto& point = iterator(point_);
  for (;;) {
    ++point;
    if (point == segment_->timeline_end()) {
      ++segment_;
      point_ = AtSegmentBegin{};
      break;
    } else if (point->first != previous_time_) {
      break;
    }
  }
  return *this;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>&
DiscreteTrajectoryIterator<Frame>::operator--() {
  NormalizeAtSegmentRBegin(point_, previous_time_);
  for (;;) {
    if (AtBegin()) {
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
  NormalizeAtSegmentTips(point, time);
  return *iterator(point);
}

template<typename Frame>
typename internal_discrete_trajectory_types::Timeline<Frame>::value_type const*
DiscreteTrajectoryIterator<Frame>::operator->() const {
  auto point = point_;
  Instant time;
  NormalizeAtSegmentTips(point, time);
  if (time == previous_time_) {
    DiscreteTrajectoryIterator it(segment_, point);
    ++it;
    point = it.point_;
  }
  return &*iterator(point);
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>::DiscreteTrajectoryIterator(
    DiscreteTrajectorySegmentIterator<Frame> const segment,
    LazyTimelineConstIterator const point)
    : segment_(segment),
      point_(point) {}

template<typename Frame>
void DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentBegin(
    LazyTimelineConstIterator& point,
    Instant& time) const {
  if (std::holds_alternative<AtSegmentBegin>(point)) {
    point = segment_->timeline_begin();
  }
  if (!std::holds_alternative<AtSegmentRBegin>(point)) {
    time = iterator(point)->first;
  }
}

template<typename Frame>
void DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentRBegin(
    LazyTimelineConstIterator& point,
    Instant& time) const {
  if (std::holds_alternative<AtSegmentRBegin>(point)) {
    point = --segment_->timeline_end();
  }
  if (!std::holds_alternative<AtSegmentBegin>(point)) {
    time = iterator(point)->first;
  }
}

template<typename Frame>
void DiscreteTrajectoryIterator<Frame>::NormalizeAtSegmentTips(
    LazyTimelineConstIterator& point,
    Instant& time) const {
  if (std::holds_alternative<AtSegmentBegin>(point)) {
    point = segment_->timeline_begin();
  }
  if (std::holds_alternative<AtSegmentRBegin>(point)) {
    point = --segment_->timeline_end();
  }
  time = iterator(point)->first;
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::Timeline::const_iterator&
DiscreteTrajectoryIterator<Frame>::iterator(
    LazyTimelineConstIterator& point) {
  return std::get<typename Timeline::const_iterator>(point);
}

template<typename Frame>
bool DiscreteTrajectoryIterator<Frame>::AtBegin() const {
  return std::holds_alternative<AtSegmentBegin>(point_) ||
         (std::holds_alternative<typename Timeline::const_iterator>(point_) &&
          iterator(point_) == segment_->timeline_begin());
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
