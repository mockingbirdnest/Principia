#pragma once

#include "testing_utilities/trajectory_factories.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_trajectory_factories {

template<typename Frame>
DiscreteTrajectory<Frame> MakeCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    Instant const& t1,
    Instant const& t2) {
  return MakeCircularTrajectory(
      ω, r, Δt, DoublePrecision<Instant>(t1), t2, trajectory);
}

template<typename Frame>
DiscreteTrajectory<Frame> AppendCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    DoublePrecision<Instant> const& t1,
    DoublePrecision<Instant> const& t2,
    DoublePrecision<Instant>& t_max) {
  Speed const v = ω * r / Radian;
  auto t = t1;
  for (; t.value <= t2; t.Increment(Δt)) {
    DegreesOfFreedom<World> const dof = {
        World::origin + Displacement<World>{{r * Cos(ω * (t.value - t0_)),
                                             r * Sin(ω * (t.value - t0_)),
                                             0 * Metre}},
        Velocity<World>{{-v * Sin(ω * (t.value - t0_)),
                         v * Cos(ω * (t.value - t0_)),
                         0 * Metre / Second}}};
    trajectory.Append(t.value, dof);
  }
  return t;
}

}  // namespace internal_trajectory_factories
}  // namespace testing_utilities
}  // namespace principia
