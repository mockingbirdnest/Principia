#pragma once

#include "testing_utilities/discrete_trajectory_factories.hpp"

#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_discrete_trajectory_factories {

using base::check_not_null;
using base::make_not_null_unique;
using geometry::Displacement;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectorySegment;
using physics::DiscreteTrajectorySegmentIterator;
using quantities::Cos;
using quantities::Sin;
using quantities::Speed;
using quantities::si::Radian;

template<typename Frame>
absl::Status DiscreteTrajectoryFactoriesFriend<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom,
    DiscreteTrajectorySegment<Frame>& segment) {
  return segment.Append(t, degrees_of_freedom);
}

template<typename Frame>
Timeline<Frame> NewMotionlessTrajectoryTimeline(Position<Frame> const& position,
                                                Time const& Δt,
                                                Instant const& t1,
                                                Instant const& t2) {
  return NewLinearTrajectoryTimeline(
      DegreesOfFreedom<Frame>(position, Velocity<Frame>()), Δt, t1, t2);
}

template<typename Frame>
Timeline<Frame>
NewLinearTrajectoryTimeline(DegreesOfFreedom<Frame> const& degrees_of_freedom,
                            Time const& Δt,
                            Instant const& t0,
                            Instant const& t1,
                            Instant const& t2) {
  Timeline<Frame> timeline;
  for (auto t = t1; t < t2; t += Δt) {
    auto const velocity = degrees_of_freedom.velocity();
    auto const position = degrees_of_freedom.position() + velocity * (t - t0);
    timeline.emplace(t, DegreesOfFreedom<Frame>(position, velocity));
  }
  return timeline;
}

template<typename Frame>
Timeline<Frame>
NewLinearTrajectoryTimeline(DegreesOfFreedom<Frame> const& degrees_of_freedom,
                            Time const& Δt,
                            Instant const& t1,
                            Instant const& t2) {
  return NewLinearTrajectoryTimeline(degrees_of_freedom, Δt, /*t0=*/t1, t1, t2);
}

template<typename Frame>
Timeline<Frame>
NewLinearTrajectoryTimeline(Velocity<Frame> const& v,
                            Time const& Δt,
                            Instant const& t1,
                            Instant const& t2) {
  return NewLinearTrajectoryTimeline(
      DegreesOfFreedom<Frame>(Frame::origin, v), Δt, /*t0=*/t1, t1, t2);
}

template<typename Frame>
Timeline<Frame>
NewCircularTrajectoryTimeline(AngularFrequency const& ω,
                              Length const& r,
                              Time const& Δt,
                              Instant const& t1,
                              Instant const& t2) {
  static Instant const t0;
  Timeline<Frame> timeline;
  Speed const v = ω * r / Radian;
  for (auto t = t1; t < t2; t += Δt) {
    DegreesOfFreedom<Frame> const dof = {
        Frame::origin + Displacement<Frame>{{r * Cos(ω * (t - t0)),
                                             r * Sin(ω * (t - t0)),
                                             Length{}}},
        Velocity<Frame>{{-v * Sin(ω * (t - t0)),
                         v * Cos(ω * (t - t0)),
                         Speed{}}}};
    timeline.emplace(t, dof);
  }
  return timeline;
}

template<typename Frame>
Timeline<Frame>
NewCircularTrajectoryTimeline(Time const& period,
                              Length const& r,
                              Time const& Δt,
                              Instant const& t1,
                              Instant const& t2) {
  return NewCircularTrajectoryTimeline<Frame>(/*ω=*/2 * π * Radian / period,
                                              r,
                                              Δt,
                                              t1,
                                              t2);
}

template<typename Frame>
void AppendTrajectoryTimeline(Timeline<Frame> const& from,
                              DiscreteTrajectorySegment<Frame>& to) {
  for (auto const& [t, degrees_of_freedom] : from) {
    EXPECT_OK(DiscreteTrajectoryFactoriesFriend<Frame>::Append(
        t, degrees_of_freedom, to));
  }
}

template<typename Frame>
void AppendTrajectoryTimeline(Timeline<Frame> const& from,
                              DiscreteTrajectory<Frame>& to) {
  for (auto const& [t, degrees_of_freedom] : from) {
    EXPECT_OK(to.Append(t, degrees_of_freedom));
  }
}

template<typename Frame>
void AppendTrajectoryTimeline(
    Timeline<Frame> const& from,
    std::function<void(
        Instant const& time,
        DegreesOfFreedom<Frame> const& degrees_of_freedom)> const& append_to) {
  for (auto const& [t, degrees_of_freedom] : from) {
    append_to(t, degrees_of_freedom);
  }
}

}  // namespace internal_discrete_trajectory_factories
}  // namespace testing_utilities
}  // namespace principia
