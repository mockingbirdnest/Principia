#pragma once

#include <list>
#include <memory>

#include "absl/container/btree_map.h"
#include "base/macros.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"

// An internal header to avoid replicating data structures in multiple places.
// Doesn't export anything outside of its internal namespace.
namespace principia {
namespace physics {

FORWARD_DECLARE_FROM(discrete_trajectory_segment,
                     TEMPLATE(typename Frame) class,
                     DiscreteTrajectorySegment);

namespace internal_discrete_trajectory_types {

using geometry::Instant;
using physics::DegreesOfFreedom;
using quantities::Length;

// |max_dense_intervals| is the maximal number of dense intervals before
// downsampling occurs.  |tolerance| is the tolerance for downsampling with
// |FitHermiteSpline|.
struct DownsamplingParameters {
  std::int64_t max_dense_intervals;
  Length tolerance;
};

// The use of an unique_ptr here makes it possible to only depend on a forward
// declaration of DiscreteTrajectorySegment.
template<typename Frame>
using Segments = std::list<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>;

template<typename Frame>
using Timeline = absl::btree_map<Instant, DegreesOfFreedom<Frame>>;

}  // namespace internal_discrete_trajectory_types
}  // namespace physics
}  // namespace principia
