#pragma once

#include <memory>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_discrete_trajectory_factories {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::DiscreteTrajectorySegment;
using physics::internal_discrete_trajectory_types::Timeline;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Time;

// A helper class to avoid having to declare friendship for all the individual
// functions in this file.
template<typename Frame>
class DiscreteTrajectoryFactoriesFriend {
 public:
  static absl::Status Append(Instant const& t,
                             DegreesOfFreedom<Frame> const& degrees_of_freedom,
                             DiscreteTrajectorySegment<Frame>& segment);
};

// A motionless trajectory at the given position.
template<typename Frame>
Timeline<Frame>
NewMotionlessTrajectoryTimeline(Position<Frame> const& position,
                                Time const& Δt,
                                Instant const& t1,
                                Instant const& t2);

// A linear trajectory with constant velocity, going through
// |degrees_of_freedom.position()| at t = t0.  The first point is at time |t1|,
// the last point at a time < |t2|.
template<typename Frame>
Timeline<Frame>
NewLinearTrajectoryTimeline(DegreesOfFreedom<Frame> const& degrees_of_freedom,
                            Time const& Δt,
                            Instant const& t0,
                            Instant const& t1,
                            Instant const& t2);
// Same as above, going through |degrees_of_freedom.position()| at t = t1.
template<typename Frame>
Timeline<Frame>
NewLinearTrajectoryTimeline(DegreesOfFreedom<Frame> const& degrees_of_freedom,
                            Time const& Δt,
                            Instant const& t1,
                            Instant const& t2);
// Same as above, going through the origin at t = t1.
template<typename Frame>
Timeline<Frame>
NewLinearTrajectoryTimeline(Velocity<Frame> const& v,
                            Time const& Δt,
                            Instant const& t1,
                            Instant const& t2);

// A trajectory with a uniform acceleration, having the specified
// |degrees_of_freedom| at t = 0.  The first point is at time |t1|, the last
// point at a time < |t2|.
template<typename Frame>
Timeline<Frame> NewAcceleratedTrajectoryTimeline(
    DegreesOfFreedom<Frame> const& degrees_of_freedom,
    Vector<Acceleration, Frame> const& acceleration,
    Time const& Δt,
    Instant const& t1,
    Instant const& t2);

// A circular trajectory in the plane XY, centred at the origin.  The first
// point is at time |t1|, the last point at a time < |t2|.
template<typename Frame>
Timeline<Frame>
NewCircularTrajectoryTimeline(AngularFrequency const& ω,
                              Length const& r,
                              Time const& Δt,
                              Instant const& t1,
                              Instant const& t2);
template<typename Frame>
Timeline<Frame>
NewCircularTrajectoryTimeline(Time const& period,
                              Length const& r,
                              Time const& Δt,
                              Instant const& t1,
                              Instant const& t2);

template<typename Frame>
void AppendTrajectoryTimeline(Timeline<Frame> const& from,
                              DiscreteTrajectorySegment<Frame>& to);

template<typename Frame>
void AppendTrajectoryTimeline(Timeline<Frame> const& from,
                              DiscreteTrajectory<Frame>& to);

template<typename Frame>
void AppendTrajectoryTimeline(
    Timeline<Frame> const& from,
    std::function<void(
        Instant const& time,
        DegreesOfFreedom<Frame> const& degrees_of_freedom)> const& append_to);

}  // namespace internal_discrete_trajectory_factories

using internal_discrete_trajectory_factories::AppendTrajectoryTimeline;
using internal_discrete_trajectory_factories::NewAcceleratedTrajectoryTimeline;
using internal_discrete_trajectory_factories::NewCircularTrajectoryTimeline;
using internal_discrete_trajectory_factories::NewLinearTrajectoryTimeline;
using internal_discrete_trajectory_factories::NewMotionlessTrajectoryTimeline;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/discrete_trajectory_factories_body.hpp"
