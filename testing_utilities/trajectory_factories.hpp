#pragma once

#include "geometry/named_quantities.hpp"
#include "numerics/double_precision.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_trajectory_factories {

using geometry::Instant;
using numerics::DoublePrecision;
using physics::DiscreteTrajectory;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Time;

template<typename Frame>
DiscreteTrajectory<Frame> MakeCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    Instant const& t1,
    Instant const& t2);

template<typename Frame>
DiscreteTrajectory<Frame> MakeCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    DoublePrecision<Instant> const& t1,
    DoublePrecision<Instant> const& t2,
    DoublePrecision<Instant>& t_max);

template<typename Frame>
void AppendTrajectory(DiscreteTrajectory<Frame> const& from,
                      DiscreteTrajectory<Frame> const& to);

}  // namespace internal_trajectory_factories
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/trajectory_factories_body.hpp"