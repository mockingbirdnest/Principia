#pragma once

#include "physics/ephemeris.hpp"

#include "base/map_util.hpp"
#include "physics/continuous_trajectory.hpp"

namespace principia {

using base::FindOrDie;

namespace physics {

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody>>> bodies,
    std::vector<DegreesOfFreedom<Frame>> initial_state,
    Instant const& initial_time,
    FixedStepSizeIntegrator<PlanetaryMotion> const& planetary_integrator,
    Time const& step_size,
    Length const& low_fitting_tolerance,
    Length const& high_fitting_tolerance)
    : bodies_(std::move(bodies)),
      planetary_integrator_(planetary_integrator),
      step_size_(step_size),
      low_fitting_tolerance_(low_fitting_tolerance),
      high_fitting_tolerance_(high_fitting_tolerance),
      last_state_(initial_state) {
  CHECK(!bodies.empty());
  CHECK_EQ(bodies.size(), initial_state.size());
  for (auto const& body : bodies_) {
    trajectories_.emplace(body.get(),
                          ContinuousTrajectory<Frame>(step_size_,
                                                      low_fitting_tolerance_,
                                                      high_fitting_tolerance_));
  }
}

template<typename Frame>
ContinuousTrajectory<Frame> const& Ephemeris<Frame>::trajectory(
    not_null<MassiveBody const*>) const {
  return FindOrDie(trajectories_, body);
}

template<typename Frame>
Instant Ephemeris<Frame>::t_min() const {
  Time t_min;
  for (auto const& pair : trajectories_) {
    ContinuousTrajectory<Frame> const& trajectory = pair.first;
    t_min = std::max(t_min, trajectory.t_min());
  }
  return t_min;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_max() const {
  Time t_max = trajectories_.begin()->first.t_max();
  for (auto const& pair : trajectories_) {
    ContinuousTrajectory<Frame> const& trajectory = pair.first;
    t_max = std::min(t_max, trajectory.t_max());
  }
  return t_max;
}

template<typename Frame>
void Ephemeris<Frame>::ForgetBefore(Instant const& t) {
  for (auto const& pair : trajectories_) {
    ContinuousTrajectory<Frame>& trajectory = pair.first;
    trajectory.ForgetBefore(t);
  }
}

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
