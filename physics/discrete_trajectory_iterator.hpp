#pragma once

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

  typename DiscreteTrajectory<Frame>::value_type const& operator*() const;
  typename DiscreteTrajectory<Frame>::value_type const* operator->() const;

 private:
  using Timeline = internal_discrete_trajectory_types::Timeline<Frame>;

  DiscreteTrajectorySegmentIterator segment_;
  Timeline::const_iterator point_;
};

}  // namespace internal_discrete_trajectory_iterator

template<typename Frame>
using DiscreteTrajectoryIterator =
    internal_discrete_trajectory_iterator::DiscreteTrajectoryIterator;

}  // namespace physics
}  // namespace principia
