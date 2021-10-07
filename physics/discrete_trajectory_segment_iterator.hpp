#pragma once

#include "absl/container/btree_map.h"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

FORWARD_DECLARE_FROM(discrete_trajectory_segment,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectorySegment);

class DiscreteTrajectoryIteratorTest;
class DiscreteTrajectorySegmentIteratorTest;

namespace internal_discrete_trajectory_segment_iterator {

using base::not_null;

template<typename Frame>
class DiscreteTrajectorySegmentIterator {
 public:
  DiscreteTrajectorySegmentIterator() = default;

  DiscreteTrajectorySegmentIterator& operator++();
  DiscreteTrajectorySegmentIterator& operator--();
  DiscreteTrajectorySegmentIterator operator++(int);
  DiscreteTrajectorySegmentIterator operator--(int);

  internal_discrete_trajectory_segment::DiscreteTrajectorySegment<Frame> const&
  operator*() const;
  internal_discrete_trajectory_segment::DiscreteTrajectorySegment<Frame> const*
  operator->() const;

 private:
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  DiscreteTrajectorySegmentIterator(not_null<Segments const*> segments,
                                    typename Segments::const_iterator iterator);
  // Not not_null<> to be default-constructible.
  Segments const* segments_ = nullptr;
  typename Segments::const_iterator iterator_;

  friend class DiscreteTrajectoryIteratorTest;
  friend class DiscreteTrajectorySegmentIteratorTest;
};

}  // namespace internal_discrete_trajectory_segment_iterator

using internal_discrete_trajectory_segment_iterator::
      DiscreteTrajectorySegmentIterator;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_segment_iterator_body.hpp"
#