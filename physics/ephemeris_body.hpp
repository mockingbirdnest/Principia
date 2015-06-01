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

  last_state_.time = initial_time;

  for (int i = 0; i < bodies_.size(); ++i) {
    auto& body = bodies_[i];
    DegreesOfFreedom<Frame> const& degrees_of_freedom = initial_state[i];
    ContinuousTrajectory<Frame>* trajectory = nullptr;
    if (body->is_oblate()) {
      // Inserting at the beginning of the vectors is O(N).
      bodies_.insert(bodies_.begin(), std::move(body));
      oblate_trajectories_.emplace(oblate_trajectories_.begin(),
                                   bodies_.front().get(),
                                   step_,
                                   low_fitting_tolerance_,
                                   high_fitting_tolerance_);
      last_state_.positions.insert(last_state_.positions.begin(),
                                   degrees_of_freedom.position());
      last_state_.velocities.insert(last_state_.velocities.begin(),
                                    degrees_of_freedom.velocity());
      trajectory = &oblate_trajectories_.front();
    } else {
      // Inserting at the end of the vectors is O(1).
      bodies_.push_back(std::move(body));
      spherical_trajectories_.emplace_back(bodies_.back().get(),
                                           step_,
                                           low_fitting_tolerance_,
                                           high_fitting_tolerance_);
      last_state_.positions.push_back(degrees_of_freedom.position());
      last_state_.velocities.push_back(degrees_of_freedom.velocity());
      trajectory = &spherical_trajectories_.back();
    }
    bodies_to_trajectories_[body.get()] = trajectory;
  }

  equation_.compute_acceleration = std::bind(somethingorother);
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
  problem.append_state = std::bind(&Ephemeris<Frame>::AppendState, this);
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
  last_state_ = state;
  int index = 0;
  for (auto const& trajectory : oblate_trajectories_) {
    trajectory.append(
        state.time,
        DegreesOfFreedom<Frame>(state.positions[index],
                                state.velocities[index]));
    ++index;
  }
  for (auto const& trajectory : spherical_trajectories_) {
    trajectory.append(
        state.time,
        DegreesOfFreedom<Frame>(state.positions[index],
                                state.velocities[index]));
    ++index;
  }
}

}  // namespace physics
}  // namespace principia
