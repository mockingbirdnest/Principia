#pragma once

#include "testing_utilities/discrete_trajectory_factories.hpp"

#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

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
DiscreteTrajectorySegment<Frame>
DiscreteTrajectoryFactoriesFriend<Frame>::MakeDiscreteTrajectorySegment(
    Segments<Frame> const& segments,
    typename Segments<Frame>::const_iterator const iterator) {
  return DiscreteTrajectorySegment<Frame>(
      DiscreteTrajectorySegmentIterator<Frame>(check_not_null(&segments),
                                               iterator));
}

template<typename Frame>
not_null<std::unique_ptr<Segments<Frame>>>
NewLinearTrajectorySegment(DegreesOfFreedom<Frame> const& degrees_of_freedom,
                           Time const& Δt,
                           Instant const& t1,
                           Instant const& t2) {
  static Instant const t0;
  auto segments = std::make_unique<Segments<Frame>>(1);
  auto const it = segments->begin();
  *it = std::make_unique<DiscreteTrajectorySegment<Frame>>(
      DiscreteTrajectoryFactoriesFriend<Frame>::MakeDiscreteTrajectorySegment(
          *segments, it));
  auto& segment = **segments->cbegin();
  for (auto t = t1; t < t2; t += Δt) {
    auto const velocity = degrees_of_freedom.velocity();
    auto const position = degrees_of_freedom.position() + velocity * (t - t0);
    DiscreteTrajectoryFactoriesFriend<Frame>::Append(
        t, DegreesOfFreedom<Frame>(position, velocity), segment);
  }
  return segments;
}

template<typename Frame>
not_null<std::unique_ptr<Segments<Frame>>>
NewLinearTrajectorySegment(Velocity<Frame> const& v,
                           Time const& Δt,
                           Instant const& t1,
                           Instant const& t2) {
  return NewLinearTrajectory(
      DegreesOfFreedom<Frame>(Frame::origin, v), Δt, t1, t2);
}

template<typename Frame>
not_null<std::unique_ptr<Segments<Frame>>>
NewCircularTrajectorySegment(AngularFrequency const& ω,
                             Length const& r,
                             Time const& Δt,
                             Instant const& t1,
                             Instant const& t2) {
  static Instant const t0;
  auto segments = std::make_unique<Segments<Frame>>(1);
  auto const it = segments->begin();
  *it = std::make_unique<DiscreteTrajectorySegment<Frame>>(
      DiscreteTrajectoryFactoriesFriend<Frame>::MakeDiscreteTrajectorySegment(
          *segments, it));
  auto& segment = **segments->cbegin();
  Speed const v = ω * r / Radian;
  for (auto t = t1; t < t2; t += Δt) {
    DegreesOfFreedom<Frame> const dof = {
        Frame::origin + Displacement<Frame>{{r * Cos(ω * (t - t0)),
                                             r * Sin(ω * (t - t0)),
                                             Length{}}},
        Velocity<Frame>{{-v * Sin(ω * (t - t0)),
                         v * Cos(ω * (t - t0)),
                         Speed{}}}};
    DiscreteTrajectoryFactoriesFriend<Frame>::Append(t, dof, segment);
  }
  return segments;
}

template<typename Frame>
not_null<std::unique_ptr<Segments<Frame>>>
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
    DiscreteTrajectoryFactoriesFriend<Frame>::Append(t, dof, to);
  }
}

}  // namespace internal_discrete_trajectory_factories
}  // namespace testing_utilities
}  // namespace principia
