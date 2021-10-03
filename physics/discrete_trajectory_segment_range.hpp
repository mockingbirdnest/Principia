#pragma once

#include <list>
#include <memory>

#include "absl/container/btree_map.h"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment_range {

using geometry::Instant;
using physics::DegreesOfFreedom;

template<typename Iterator>
class DiscreteTrajectorySegmentRange {
 public:
  DiscreteTrajectorySegmentRange() = default;

  Iterator begin() const;
  Iterator end() const;

 private:
  Iterator begin_;
  Iterator end_;
};

}  // namespace internal_discrete_trajectory_segment_range

template<typename Iterator>
using DiscreteTrajectorySegmentRange =
    internal_discrete_trajectory_segment_range::DiscreteTrajectorySegmentRange;

}  // namespace physics
}  // namespace principia
