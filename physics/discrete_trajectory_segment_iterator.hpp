#pragma once

#include <iterator>

#include "absl/container/btree_map.h"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "physics/discrete_trajectory_segment_range.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {

namespace testing_utilities {
FORWARD_DECLARE_FROM(discrete_trajectory_factories,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectoryFactoriesFriend);
}  // namespace testing_utilities

namespace physics {

FORWARD_DECLARE_FROM(discrete_trajectory,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectory);
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
  using value_type = DiscreteTrajectorySegment<Frame>;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::bidirectional_iterator_tag;

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

  DiscreteTrajectorySegmentIterator(not_null<Segments*> segments,
                                    typename Segments::iterator iterator);

  bool is_begin() const;
  bool is_end() const;
  DiscreteTrajectorySegmentRange<DiscreteTrajectorySegmentIterator>
  segments() const;

  typename Segments::iterator iterator() const;

  // Not not_null<> to be default-constructible.
  Segments* segments_ = nullptr;
  typename Segments::iterator iterator_;

  template<typename F>
  friend class physics::DiscreteTrajectory;
  template<typename F>
  friend class physics::DiscreteTrajectoryIterator;

  // For testing.
  friend class physics::DiscreteTrajectoryIteratorTest;
  friend class physics::DiscreteTrajectorySegmentIteratorTest;
  friend class physics::DiscreteTrajectorySegmentTest;
  template<typename F>
  friend class testing_utilities::DiscreteTrajectoryFactoriesFriend;
};

}  // namespace internal_discrete_trajectory_segment_iterator

using internal_discrete_trajectory_segment_iterator::
      DiscreteTrajectorySegmentIterator;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_segment_iterator_body.hpp"
