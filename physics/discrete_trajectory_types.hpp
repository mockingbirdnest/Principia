#pragma once

#include <list>

#include "absl/container/btree_set.h"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "geometry/instant.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"

// An internal header to avoid replicating data structures in multiple places.
// Doesn't export anything outside of its internal namespace.
namespace principia {
namespace physics {

FORWARD_DECLARE(TEMPLATE(typename Frame) class,
                DiscreteTrajectorySegment,
                FROM(discrete_trajectory_segment),
                INTO(discrete_trajectory_types));

namespace _discrete_trajectory_types {
namespace internal {

using namespace principia::geometry::_instant;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_quantities;

// `max_dense_intervals` is the maximal number of dense intervals before
// downsampling occurs.  `tolerance` is the tolerance for downsampling with
// `FitHermiteSpline`.
struct DownsamplingParameters {
  std::int64_t max_dense_intervals;
  Length tolerance;
};

template<typename Frame>
struct value_type {
  value_type(Instant const& time,
             DegreesOfFreedom<Frame> const& degrees_of_freedom);
  Instant time;
  DegreesOfFreedom<Frame> degrees_of_freedom;
};

struct Earlier {
  using is_transparent = void;

  Earlier() = default;

  template<typename Frame>
  bool operator()(value_type<Frame> const& left,
                  value_type<Frame> const& right) const;
  template<typename Frame>
  bool operator()(Instant const& left, value_type<Frame> const& right) const;
  template<typename Frame>
  bool operator()(value_type<Frame> const& left, Instant const& right) const;
};

template<typename Frame>
using Segments = std::list<DiscreteTrajectorySegment<Frame>>;

template<typename Frame>
using Timeline = absl::btree_set<value_type<Frame>, Earlier>;

}  // namespace internal

using internal::DownsamplingParameters;
using internal::Segments;
using internal::Timeline;

}  // namespace _discrete_trajectory_types
}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_types_body.hpp"
