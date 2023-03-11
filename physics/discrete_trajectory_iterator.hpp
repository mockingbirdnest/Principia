#pragma once

#include <cstdint>
#include <iterator>
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

namespace _discrete_trajectory_iterator {
namespace internal {

using namespace principia::geometry::_named_quantities;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory_segment;
using namespace principia::physics::_discrete_trajectory_types;

template<typename Frame>
class DiscreteTrajectoryIterator {
 public:
  using difference_type = std::int64_t;
  using value_type = typename Timeline<Frame>::value_type;
  using pointer = value_type const*;
  using reference = value_type const&;
  using iterator_category = std::random_access_iterator_tag;

  DiscreteTrajectoryIterator() = default;

  DiscreteTrajectoryIterator& operator++();
  DiscreteTrajectoryIterator& operator--();
  DiscreteTrajectoryIterator operator++(int);
  DiscreteTrajectoryIterator operator--(int);

  reference operator*() const;
  pointer operator->() const;

  DiscreteTrajectoryIterator& operator+=(difference_type n);
  DiscreteTrajectoryIterator& operator-=(difference_type n);
  reference operator[](difference_type n) const;

  // The operator+ are outside of this class because one of them cannot be a
  // member.
  DiscreteTrajectoryIterator operator-(difference_type n) const;
  difference_type operator-(DiscreteTrajectoryIterator right) const;

  bool operator==(DiscreteTrajectoryIterator other) const;
  bool operator!=(DiscreteTrajectoryIterator other) const;
  bool operator<(DiscreteTrajectoryIterator other) const;
  bool operator>(DiscreteTrajectoryIterator other) const;
  bool operator<=(DiscreteTrajectoryIterator other) const;
  bool operator>=(DiscreteTrajectoryIterator other) const;

 private:
  using Timeline = _discrete_trajectory_types::Timeline<Frame>;

  // Optional because we cannot construct a point iterator in the end segment.
  using OptionalTimelineConstIterator =
      std::optional<typename Timeline::const_iterator>;

  // Constructs an `end()` iterator.
  static DiscreteTrajectoryIterator EndOfLastSegment(
      DiscreteTrajectorySegmentIterator<Frame> segment);

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
  friend class _discrete_trajectory::internal::DiscreteTrajectory;
  template<typename F>
  friend class _discrete_trajectory_segment::internal::
      DiscreteTrajectorySegment;
};

template<typename Frame>
DiscreteTrajectoryIterator<Frame> operator+(
    DiscreteTrajectoryIterator<Frame> it,
    typename DiscreteTrajectoryIterator<Frame>::difference_type n);
template<typename Frame>
DiscreteTrajectoryIterator<Frame> operator+(
    typename DiscreteTrajectoryIterator<Frame>::difference_type n,
    DiscreteTrajectoryIterator<Frame> it);

}  // namespace internal

using internal::DiscreteTrajectoryIterator;

}  // namespace _discrete_trajectory_iterator
}  // namespace physics
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_discrete_trajectory_iterator;
}  // namespace principia::physics

#include "physics/discrete_trajectory_iterator_body.hpp"
