#pragma once

#include <iterator>

#include "absl/container/btree_map.h"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {

namespace testing_utilities {
FORWARD_DECLARE(TEMPLATE(typename Frame) class,
                DiscreteTrajectoryFactoriesFriend,
                FROM(discrete_trajectory_factories));
}  // namespace testing_utilities

namespace physics {

FORWARD_DECLARE(TEMPLATE(typename Frame) class,
                DiscreteTrajectory,
                FROM(discrete_trajectory),
                INTO(discrete_trajectory_segment_iterator));
FORWARD_DECLARE(TEMPLATE(typename Frame) class,
                DiscreteTrajectoryIterator,
                FROM(discrete_trajectory_iterator),
                INTO(discrete_trajectory_segment_iterator));
FORWARD_DECLARE(TEMPLATE(typename Frame) class,
                DiscreteTrajectorySegment,
                FROM(discrete_trajectory_segment),
                INTO(discrete_trajectory_segment_iterator));

class DiscreteTrajectoryIteratorTest;
class DiscreteTrajectorySegmentIteratorTest;
class DiscreteTrajectorySegmentTest;

namespace _discrete_trajectory_segment_iterator {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::physics::_discrete_trajectory_types;

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
  DiscreteTrajectorySegmentIterator EndSegment() const;

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

#include "physics/discrete_trajectory_segment_iterator_body.hpp"
