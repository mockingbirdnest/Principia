#pragma once

#include <optional>

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

  // Optional because we cannot construct a point iterator in the end segment.
  using OptionalTimelineConstIterator =
      std::optional<typename Timeline::const_iterator>;

  DiscreteTrajectoryIterator(DiscreteTrajectorySegmentIterator<Frame> segment,
                             OptionalTimelineConstIterator point);

  static bool is_at_end(OptionalTimelineConstIterator point);

  static typename Timeline::const_iterator& iterator(
      OptionalTimelineConstIterator& point);
  static typename Timeline::const_iterator const& iterator(
      OptionalTimelineConstIterator const& point);

  // |point_| is always an iterator in the timeline of the segment denoted by
  // |segment_|.  When |segment_| is at the end of its list, |point_| is
  // nullopt.  It is possible to have repeated times in a segment or across
  // segments and the iterator will skip them, so that they will appear as a
  // single point to clients.
  DiscreteTrajectorySegmentIterator<Frame> segment_;
  OptionalTimelineConstIterator point_;

  // The last time that was seen by the iterator.  Used to skip over repeated
  // times.
  std::optional<Instant> previous_time_;

  friend class DiscreteTrajectoryIteratorTest;
};

}  // namespace internal_discrete_trajectory_iterator

using internal_discrete_trajectory_iterator::DiscreteTrajectoryIterator;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_iterator_body.hpp"
