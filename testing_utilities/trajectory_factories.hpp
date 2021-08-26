#pragma once

#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/double_precision.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_trajectory_factories {

using base::not_null;
using geometry::Instant;
using numerics::DoublePrecision;
using physics::DiscreteTrajectory;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Time;

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>> MakeCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    Instant const& t1,
    Instant const& t2);

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>> MakeCircularTrajectory(
    AngularFrequency const& ω,
    Length const& r,
    Time const& Δt,
    DoublePrecision<Instant> const& t1,
    DoublePrecision<Instant> const& t2,
    DoublePrecision<Instant>& t_max);

template<typename Frame>
void AppendTrajectory(DiscreteTrajectory<Frame> const& from,
                      DiscreteTrajectory<Frame>& to);

}  // namespace internal_trajectory_factories

using internal_trajectory_factories::AppendTrajectory;
using internal_trajectory_factories::MakeCircularTrajectory;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/trajectory_factories_body.hpp"
