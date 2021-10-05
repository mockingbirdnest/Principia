#pragma once

#include <optional>
#include <variant>

#include "absl/container/btree_map.h"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

template<typename Frame>
class DiscreteTrajectory;

namespace internal_discrete_trajectory_iterator {

using geometry::Instant;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectoryIterator {
 public:
  DiscreteTrajectoryIterator() = default;

  DiscreteTrajectoryIterator& operator++();
  DiscreteTrajectoryIterator& operator--();
  DiscreteTrajectoryIterator operator++(int);
  DiscreteTrajectoryIterator operator--(int);

  typename
  internal_discrete_trajectory_types::Timeline<Frame>::value_type const&
  operator*() const;
  typename
  internal_discrete_trajectory_types::Timeline<Frame>::value_type const*
  operator->() const;

  bool operator==(DiscreteTrajectoryIterator const& other) const;
  bool operator!=(DiscreteTrajectoryIterator const& other) const;

 private:
  using Timeline = internal_discrete_trajectory_types::Timeline<Frame>;

  //TODO(phl):comment
  struct AtSegmentBegin {};
  struct AtSegmentRBegin {};
  using LazyTimelineConstIterator =
      std::variant<AtSegmentBegin,
                   AtSegmentRBegin,
                   typename Timeline::const_iterator>;

  DiscreteTrajectoryIterator(DiscreteTrajectorySegmentIterator<Frame> segment,
                             LazyTimelineConstIterator point);

  bool IsAtSegmentBegin() const;
  bool IsAtSegmentRBegin() const;

  static void NormalizeAtSegmentBegin(
      DiscreteTrajectorySegmentIterator<Frame> const& segment,
      LazyTimelineConstIterator& point,
      Instant& time);
  static void NormalizeAtSegmentRBegin(
      DiscreteTrajectorySegmentIterator<Frame> const& segment,
      LazyTimelineConstIterator& point,
      Instant& time);
  static void NormalizeAtSegmentTips(
      DiscreteTrajectorySegmentIterator<Frame> const& segment,
      LazyTimelineConstIterator& point,
      Instant& time);

  static typename Timeline::const_iterator& iterator(
      LazyTimelineConstIterator& point);
  static typename Timeline::const_iterator const& iterator(
      LazyTimelineConstIterator const& point);

  // |point_| is always an iterator in the timeline of the segment denoted by
  // |segment_|.
  // TODO(phl): Figure out what to do with empty segments.
  DiscreteTrajectorySegmentIterator<Frame> segment_;
  LazyTimelineConstIterator point_;

  //TODO(phl):optional?
  Instant previous_time_;

  friend class DiscreteTrajectoryIteratorTest;
};

}  // namespace internal_discrete_trajectory_iterator

using internal_discrete_trajectory_iterator::DiscreteTrajectoryIterator;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_iterator_body.hpp"
