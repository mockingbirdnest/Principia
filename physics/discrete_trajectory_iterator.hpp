#pragma once

#include <cstdint>
#include <optional>

#include "absl/container/btree_map.h"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

FORWARD_DECLARE_FROM(discrete_trajectory_segment,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectorySegment);

namespace internal_discrete_trajectory_iterator {

using geometry::Instant;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectoryIterator {
 public:
  using difference_type = std::int64_t;
  using value_type =
      typename internal_discrete_trajectory_types::Timeline<Frame>::value_type;
  using pointer = value_type const*;
  using reference = value_type const&;

  DiscreteTrajectoryIterator() = default;

  DiscreteTrajectoryIterator& operator++();
  DiscreteTrajectoryIterator& operator--();
  DiscreteTrajectoryIterator operator++(int);
  DiscreteTrajectoryIterator operator--(int);

  reference operator*() const;
  pointer operator->() const;

  DiscreteTrajectoryIterator& operator+=(difference_type n);
  DiscreteTrajectoryIterator& operator-=(difference_type n);
  friend DiscreteTrajectoryIterator operator+(DiscreteTrajectoryIterator it,
                                              difference_type n);
  friend DiscreteTrajectoryIterator operator+(difference_type n,
                                              DiscreteTrajectoryIterator it);
  friend DiscreteTrajectoryIterator operator-(DiscreteTrajectoryIterator it,
                                              difference_type n);
  reference operator[](difference_type n) const;

  bool operator==(DiscreteTrajectoryIterator other) const;
  bool operator!=(DiscreteTrajectoryIterator other) const;
  bool operator<(DiscreteTrajectoryIterator other) const;
  bool operator>(DiscreteTrajectoryIterator other) const;
  bool operator<=(DiscreteTrajectoryIterator other) const;
  bool operator>=(DiscreteTrajectoryIterator other) const;

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

  template<typename F>
  friend class physics::DiscreteTrajectorySegment;
  friend class physics::DiscreteTrajectoryIteratorTest;
};

}  // namespace internal_discrete_trajectory_iterator

using internal_discrete_trajectory_iterator::DiscreteTrajectoryIterator;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_iterator_body.hpp"
