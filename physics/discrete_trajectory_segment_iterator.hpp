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

namespace _discrete_trajectory_segment_iterator {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_iterator;
using namespace principia::physics::_discrete_trajectory_segment;

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
  using Segments = _discrete_trajectory_types::Segments<Frame>;

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
  friend class _discrete_trajectory::internal::DiscreteTrajectory;
  template<typename F>
  friend class _discrete_trajectory_iterator::internal::
      DiscreteTrajectoryIterator;

  // For testing.
  friend class physics::DiscreteTrajectoryIteratorTest;
  friend class physics::DiscreteTrajectorySegmentIteratorTest;
  friend class physics::DiscreteTrajectorySegmentTest;
  template<typename F>
  friend class testing_utilities::_discrete_trajectory_factories::
      DiscreteTrajectoryFactoriesFriend;
};

}  // namespace internal

using internal::DiscreteTrajectorySegmentIterator;

}  // namespace _discrete_trajectory_segment_iterator
}  // namespace physics
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_discrete_trajectory_segment_iterator;
}  // namespace principia::physics

#include "physics/discrete_trajectory_segment_iterator_body.hpp"
