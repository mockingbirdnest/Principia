#pragma once

#include <list>
#include <memory>

#include "absl/container/btree_map.h"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

template<typename Frame>
class DiscreteTrajectorySegment;

namespace internal_discrete_trajectory_segment_iterator {

using geometry::Instant;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectorySegmentIterator {
 public:
  DiscreteTrajectorySegmentIterator() = default;

  DiscreteTrajectorySegmentIterator& operator++();
  DiscreteTrajectorySegmentIterator& operator--();
  DiscreteTrajectorySegmentIterator operator++(int);
  DiscreteTrajectorySegmentIterator operator--(int);

  DiscreteTrajectorySegment<Frame> const& operator*() const;
  DiscreteTrajectorySegment<Frame> const* operator->() const;

 private:
  using Segments = internal_discrete_trajectory_types::Segments<Frame>;

  Segments::const_iterator segment_;
};

}  // namespace internal_discrete_trajectory_segment_iterator

template<typename Frame>
using DiscreteTrajectorySegmentIterator =
    internal_discrete_trajectory_segment_iterator::
        DiscreteTrajectorySegmentIterator;

}  // namespace principia
}  // namespace physics
