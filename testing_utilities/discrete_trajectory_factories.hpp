#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_discrete_trajectory_factories {

using base::not_null;
using geometry::Instant;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectorySegment;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Time;

// A linear trajectory with constant velocity, going through
// |degrees_of_freedom.position()| at t = 0.  The first point is at time |t1|,
// the last point at a time < |t2|.
template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewLinearTrajectorySegment(DegreesOfFreedom<Frame> const& degrees_of_freedom,
                           Time const& Δt,
                           Instant const& t1,
                           Instant const& t2);
// Same as above, going through the origin at t = 0.
template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewLinearTrajectorySegment(Velocity<Frame> const& v,
                           Time const& Δt,
                           Instant const& t1,
                           Instant const& t2);

// A circular trajectory in the plane XY, centred at the origin.  The first
// point is at time |t1|, the last point at a time < |t2|.
template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewCircularTrajectorySegment(AngularFrequency const& ω,
                             Length const& r,
                             Time const& Δt,
                             Instant const& t1,
                             Instant const& t2);
template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewCircularTrajectorySegment(Time const& period,
                             Length const& r,
                             Time const& Δt,
                             Instant const& t1,
                             Instant const& t2);

template<typename Frame>
void AppendTrajectorySegment(DiscreteTrajectorySegment<Frame> const& from,
                             DiscreteTrajectorySegment<Frame>& to);

}  // namespace internal_discrete_trajectory_factories

using internal_discrete_trajectory_factories::AppendTrajectorySegment;
using internal_discrete_trajectory_factories::NewCircularTrajectorySegment;
using internal_discrete_trajectory_factories::NewLinearTrajectorySegment;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/discrete_trajectory_factories_body.hpp"
