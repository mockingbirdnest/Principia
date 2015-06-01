#pragma once

#include "physics/ephemeris.hpp"

#include "base/map_util.hpp"
#include "physics/continuous_trajectory.hpp"

namespace principia {

using base::FindOrDie;
using integrators::IntegrationProblem;

namespace physics {

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody>>> bodies,
    std::vector<DegreesOfFreedom<Frame>> initial_state,
    Instant const& initial_time,
    FixedStepSizeIntegrator<NewtonianMotionEquation> const&
        planetary_integrator,
    Time const& step,
    Length const& low_fitting_tolerance,
    Length const& high_fitting_tolerance)
    : planetary_integrator_(planetary_integrator),
      step_(step),
      low_fitting_tolerance_(low_fitting_tolerance),
      high_fitting_tolerance_(high_fitting_tolerance) {
  CHECK(!bodies.empty());
  CHECK_EQ(bodies.size(), initial_state.size());
  for (auto& body : bodies) {
    bodies_and_trajectories_.emplace_back(
        std::move(body),
        ContinuousTrajectory<Frame>(step_,
                                    low_fitting_tolerance_,
                                    high_fitting_tolerance_));
    bodies_to_trajectories_[bodies_and_trajectories_.back().first.get()] =
        &bodies_and_trajectories_.back().second;
  }

  equation_.compute_acceleration = std::bind(somethingorother);

  last_state_.time = initial_time;
  for (auto const& degrees_of_freedom : initial_state) {
    last_state_.positions.push_back(degrees_of_freedom.position());
    last_state_.velocities.push_back(degrees_of_freedom.velocity());
  }
}

template<typename Frame>
ContinuousTrajectory<Frame> const& Ephemeris<Frame>::trajectory(
    not_null<MassiveBody const*>) const {
  return FindOrDie(bodies_to_trajectories_, body);
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
void Ephemeris<Frame>::Prolong(Instant const& t) {
  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation = equation_;
  problem.append_state = std::bind();
  problem.t_final = t;
  problem.initial_state = &last_state_;

  planetary_integrator_.Solve(problem, step_);

  //TODO(phl): Ensure that t_max has reached t.
}

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

template<typename Frame>
void Ephemeris<Frame>::AppendState(
    typename NewtonianMotionEquation::SystemState const& state) {
}

}  // namespace physics
}  // namespace principia
