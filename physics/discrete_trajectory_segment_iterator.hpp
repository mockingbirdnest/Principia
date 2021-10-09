#pragma once

#include "absl/container/btree_map.h"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

FORWARD_DECLARE_FROM(discrete_trajectory_iterator,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectoryIterator);
FORWARD_DECLARE_FROM(discrete_trajectory_segment,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectorySegment);

class DiscreteTrajectoryIteratorTest;
class DiscreteTrajectorySegmentIteratorTest;
class DiscreteTrajectorySegmentTest;

namespace internal_discrete_trajectory_segment_iterator {

using base::not_null;

template<typename Frame>
class DiscreteTrajectorySegmentIterator {
 public:
  using difference_type = std::int64_t;
  using pointer = value_type const*;
  using reference = value_type const&;
  using value_type = DiscreteTrajectorySegment<Frame>;

  DiscreteTrajectorySegmentIterator() = default;

  DiscreteTrajectorySegmentIterator& operator++();
  DiscreteTrajectorySegmentIterator& operator--();
  DiscreteTrajectorySegmentIterator operator++(int);
  DiscreteTrajectorySegmentIterator operator--(int);

  reference operator*() const;
  pointer operator->() const;

  bool operator==(DiscreteTrajectorySegmentIterator const& other) const;
  bool operator!=(DiscreteTrajectorySegmentIterator const& other) const;

 private:
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  DiscreteTrajectorySegmentIterator(not_null<Segments const*> segments,
                                    typename Segments::const_iterator iterator);

  DiscreteTrajectorySegmentIterator begin() const;
  DiscreteTrajectorySegmentIterator end() const;

  // Not not_null<> to be default-constructible.
  Segments const* segments_ = nullptr;
  typename Segments::const_iterator iterator_;

  template<typename F>
  friend class DiscreteTrajectoryIterator;
  friend class DiscreteTrajectoryIteratorTest;
  friend class DiscreteTrajectorySegmentIteratorTest;
  friend class DiscreteTrajectorySegmentTest;
};

}  // namespace internal_discrete_trajectory_segment_iterator

using internal_discrete_trajectory_segment_iterator::
      DiscreteTrajectorySegmentIterator;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_segment_iterator_body.hpp"
