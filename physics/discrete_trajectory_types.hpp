#pragma once

#include <list>
#include <memory>
#include <tuple>

#include "absl/container/btree_set.h"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "numerics/hermite3.hpp"
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
using namespace principia::geometry::_space;
using namespace principia::numerics::_hermite3;
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

  // This interpolation is not for use by clients (hence the structured bindings
  // definitions below) but is exclusively for use by
  // `DiscreteTrajectorySegment`.  It is appropriate for evalutions on the
  // left-open, right-closed trajectory interval that ends at `time`.  It is
  // missing (set to `nullptr`) for the first point of a segment.
//TODO(phl)comment
  std::shared_ptr<Hermite3<Position<Frame>, Instant>> interpolation;

  // Support for structured bindings of `time` and `degrees_of_freedom`.
  template<std::size_t i, typename Self>
  constexpr auto&& get(this Self&& self);
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

namespace std {

template<typename Frame>
struct tuple_size<
    principia::physics::_discrete_trajectory_types::internal::value_type<Frame>>
    : std::integral_constant<std::size_t, 2> {};

template<std::size_t i, typename Frame>
struct tuple_element<
    i,
    principia::physics::_discrete_trajectory_types::internal::value_type<Frame>>
    : tuple_element<i,
                    tuple<principia::geometry::_instant::Instant,
                          principia::physics::_degrees_of_freedom::
                              DegreesOfFreedom<Frame>>> {};

}  // namespace std

#include "physics/discrete_trajectory_types_body.hpp"
