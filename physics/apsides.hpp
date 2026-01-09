#pragma once

#include <functional>
#include <vector>

#include "absl/container/btree_set.h"
#include "absl/status/status.h"
#include "base/constant_function.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/rotating_body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace _apsides {
namespace internal {

using namespace principia::base::_constant_function;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_trajectory;
using namespace principia::quantities::_quantities;

template<typename Frame>
using DistinguishedPoints = absl::btree_map<Instant, DegreesOfFreedom<Frame>>;

// Computes the apsides with respect to `reference` for the section given by
// `begin` and `end` of `trajectory`.  Appends to the given output trajectories
// one point for each apsis.
template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    Trajectory<Frame> const& trajectory,
                    typename DiscreteTrajectory<Frame>::iterator begin,
                    typename DiscreteTrajectory<Frame>::iterator end,
                    Instant const& t_max,
                    int max_points,
                    DistinguishedPoints<Frame>& apoapsides,
                    DistinguishedPoints<Frame>& periapsides);

// Returns the ordered time intervals where there can be a collision with the
// `reference_body` because the `trajectory` is below its `max_radius`.
template<typename Frame>
std::vector<Interval<Instant>> ComputeCollisionIntervals(
    RotatingBody<Frame> const& reference_body,
    Trajectory<Frame> const& reference,
    Trajectory<Frame> const& trajectory,
    DistinguishedPoints<Frame> const& apoapsides,
    DistinguishedPoints<Frame> const& periapsides);

// Computes the first collision between a vessel and a rotating body over the
// given time `interval` with an accuracy better than `max_collision_error`.
// Returns `nullopt` if there is no collision over the `interval`.  The
// `interval` should have been obtained by `ComputeCollisionIntervals`. `radius`
// must give the radius of the celestial at a particular position given by its
// latitude and longitude.  It must never exceed the `max_radius` of the body.
template<typename Frame>
std::optional<typename DistinguishedPoints<Frame>::value_type>
ComputeFirstCollision(
    RotatingBody<Frame> const& reference_body,
    Trajectory<Frame> const& reference,
    Trajectory<Frame> const& trajectory,
    Interval<Instant> const& interval,
    Length const& max_collision_error,
    std::function<Length(Angle const& latitude, Angle const& longitude)> const&
        radius);

// Computes the crossings of the section given by `begin` and `end` of
// `trajectory` with the xy plane.  Appends the crossings that go towards the
// `north` side of the xy plane to `ascending`, and those that go away from the
// `north` side to `descending`.
// Nodes for which `predicate` returns false are excluded.
template<typename Frame, typename Predicate = ConstantFunction<bool>>
absl::Status ComputeNodes(Trajectory<Frame> const& trajectory,
                          typename DiscreteTrajectory<Frame>::iterator begin,
                          typename DiscreteTrajectory<Frame>::iterator end,
                          Instant const& t_max,
                          Vector<double, Frame> const& north,
                          int max_points,
                          DistinguishedPoints<Frame>& ascending,
                          DistinguishedPoints<Frame>& descending,
                          Predicate predicate = Identically(true));

// TODO(egg): when we can usefully iterate over an arbitrary `Trajectory`, move
// the following from `Ephemeris`.
#if 0
template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& trajectory1,
                    Trajectory<Frame> const& trajectory2,
                    DistinguishedPoints<Frame>& apoapsides1,
                    DistinguishedPoints<Frame>& periapsides1,
                    DistinguishedPoints<Frame>& apoapsides2,
                    DistinguishedPoints<Frame>& periapsides2);
#endif

}  // namespace internal

using internal::ComputeApsides;
using internal::ComputeCollisionIntervals;
using internal::ComputeFirstCollision;
using internal::ComputeNodes;
using internal::DistinguishedPoints;

}  // namespace _apsides
}  // namespace physics
}  // namespace principia

#include "physics/apsides_body.hpp"
