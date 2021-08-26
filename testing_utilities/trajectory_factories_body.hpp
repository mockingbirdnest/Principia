#pragma once

#include "testing_utilities/trajectory_factories.hpp"

#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_trajectory_factories {

using base::make_not_null_unique;
using geometry::Displacement;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using quantities::Cos;
using quantities::Sin;
using quantities::Speed;
using quantities::si::Radian;

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>> NewCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    Instant const& t1,
    Instant const& t2) {
  DoublePrecision<Instant> t_max;
  return NewCircularTrajectory<Frame>(ω,
                                      r,
                                      Δt,
                                      DoublePrecision<Instant>(t1),
                                      DoublePrecision<Instant>(t2),
                                      t_max);
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>> NewCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    DoublePrecision<Instant> const& t1,
    DoublePrecision<Instant> const& t2,
    DoublePrecision<Instant>& t_max) {
  static Instant const t0;
  Speed const v = ω * r / Radian;
  auto trajectory = make_not_null_unique<DiscreteTrajectory<Frame>>();
  auto t = t1;
  for (; t.value <= t2.value; t.Increment(Δt)) {
    DegreesOfFreedom<Frame> const dof = {
        Frame::origin + Displacement<Frame>{{r * Cos(ω * (t.value - t0)),
                                             r * Sin(ω * (t.value - t0)),
                                             Length{}}},
        Velocity<Frame>{{-v * Sin(ω * (t.value - t0)),
                         v * Cos(ω * (t.value - t0)),
                         Speed{}}}};
    trajectory->Append(t.value, dof);
  }
  t_max = t;
  return std::move(trajectory);
}

template<typename Frame>
void AppendTrajectory(DiscreteTrajectory<Frame> const& from,
                      DiscreteTrajectory<Frame>& to) {
  for (auto const& [t, dof] : from) {
    to.Append(t, dof);
  }
}

}  // namespace internal_trajectory_factories
}  // namespace testing_utilities
}  // namespace principia
