#pragma once

#include "physics/ephemeris.hpp"

namespace principia {
namespace physics {

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody>>> bodies,
    std::vector<DegreesOfFreedom<Frame>> initial_state,
    Instant const& initial_time,
    FixedStepSizeIntegrator<PlanetaryMotion> const& planetary_integrator,
    Time const& step_size,
    Length const& low_fitting_tolerance,
    Length const& high_fitting_tolerance) {}

template<typename Frame>
ContinuousTrajectory<Frame> const& Ephemeris<Frame>::trajectory(
    not_null<MassiveBody const*>) const {}

template<typename Frame>
Instant Ephemeris<Frame>::t_min() const {}

template<typename Frame>
Instant Ephemeris<Frame>::t_max() const {}

template<typename Frame>
void Ephemeris<Frame>::ForgetBefore(Instant const& t) {}

template<typename Frame>
void Ephemeris<Frame>::Prolong(Instant const& t) {}

template<typename Frame>
void Ephemeris<Frame>::Flow(
    not_null<Trajectory<Frame>*> const trajectory,
    std::function<
        Vector<Acceleration, Frame>(
            Instant const&)> intrinsic_acceleration,
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance,
    AdaptiveStepSizeIntegrator<TimedBurnMotion> integrator,
    Instant const& t) {}


}  // namespace physics
}  // namespace principia
