#pragma once

#include "testing_utilities/discrete_trajectory_factories.hpp"

#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_discrete_trajectory_factories {

using base::make_not_null_unique;
using geometry::Displacement;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using quantities::Cos;
using quantities::Sin;
using quantities::Speed;
using quantities::si::Radian;

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewLinearTrajectorySegment(DegreesOfFreedom<Frame> const& degrees_of_freedom,
                           Time const& Δt,
                           Instant const& t1,
                           Instant const& t2) {
  static Instant const t0;
  auto segment = make_not_null_unique<DiscreteTrajectorySegment<Frame>>();
  for (auto t = t1; t < t2; t += Δt) {
    auto const velocity = degrees_of_freedom.velocity();
    auto const position = degrees_of_freedom.position() + velocity * (t - t0);
    segment->Append(t, DegreesOfFreedom<Frame>(position, velocity));
  }
  return segment;
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewLinearTrajectorySegment(Velocity<Frame> const& v,
                           Time const& Δt,
                           Instant const& t1,
                           Instant const& t2) {
  return NewLinearTrajectory(
      DegreesOfFreedom<Frame>(Frame::origin, v), Δt, t1, t2);
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewCircularTrajectorySegment(AngularFrequency const& ω,
                             Length const& r,
                             Time const& Δt,
                             Instant const& t1,
                             Instant const& t2) {
  static Instant const t0;
  Speed const v = ω * r / Radian;
  auto segment = make_not_null_unique<DiscreteTrajectorySegment<Frame>>();
  for (auto t = t1; t < t2; t += Δt) {
    DegreesOfFreedom<Frame> const dof = {
        Frame::origin + Displacement<Frame>{{r * Cos(ω * (t - t0)),
                                             r * Sin(ω * (t - t0)),
                                             Length{}}},
        Velocity<Frame>{{-v * Sin(ω * (t - t0)),
                         v * Cos(ω * (t - t0)),
                         Speed{}}}};
    segment->Append(t, dof);
  }
  return segment;
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectorySegment<Frame>>>
NewCircularTrajectorySegment(Time const& period,
                             Length const& r,
                             Time const& Δt,
                             Instant const& t1,
                             Instant const& t2) {
  return NewCircularTrajectory<Frame>(/*ω=*/2 * π * Radian / period,
                                      r,
                                      Δt,
                                      t1,
                                      t2);
}

template<typename Frame>
void AppendTrajectorySegment(DiscreteTrajectorySegment<Frame> const& from,
                             DiscreteTrajectorySegment<Frame>& to) {
  for (auto const& [t, dof] : from) {
    to.Append(t, dof);
  }
}

}  // namespace internal_discrete_trajectory_factories
}  // namespace testing_utilities
}  // namespace principia
