#pragma once

#include <functional>
#include <map>
#include <memory>

#include "base/not_null.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "physics/trajectory.hpp"

namespace principia {

using integrators::AdaptiveStepSizeIntegrator;
using integrators::FixedStepSizeIntegrator;
using integrators::SpecialSecondOrderDifferentialEquation;

namespace physics {

template<typename Frame>
class Ephemeris {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  // The equation describing the motion of the |bodies_|.
  using PlanetaryMotion =
      SpecialSecondOrderDifferentialEquation<Position<Frame>>;
  PlanetaryMotion::SystemState last_state_;
  // We don't specify non-autonomy in PlanetaryMotion since there isn't a type
  // for that at this time, so time-dependent intrinsic acceleration yields the
  // same type of map.
  using TimedBurnMotion = PlanetaryMotion;

  // Constructs an Ephemeris that owns the |bodies|.
  Ephemeris(
      std::vector<not_null<std::unique_ptr<MassiveBody>>> bodies,
      std::vector<DegreesOfFreedom<Frame>> initial_state,
      Instant const& initial_time,
      FixedStepSizeIntegrator<PlanetaryMotion> const& planetary_integrator,
      Time const& step_size);

  ContinuousTrajectory<Frame> const& trajectory(
      not_null<MassiveBody const*>) const;

  // The common |t_min| of the trajectories.
  Instant t_min() const;
  // The common |t_max| of the trajectories.
  Instant t_max() const;

  // Calls |ForgetBefore| on all trajectories.
  void ForgetBefore(Instant const& t);

  // Prolongs the ephemeris up to at least |t|.  After the call, |t_max() >= t|.
  void Prolong(Instant const& t);

  // Computes the trajectory of a massless body in the gravitational potential
  // described by |this| and subject to the given |intrinsic_acceleration| up to
  // exactly |t|.
  // If |t > t_max()|, calls |Prolong(t)| beforehand.
  void Advance(
      not_null<Trajectory<Frame>*> const trajectory,
      std::function<
          Vector<Acceleration, Frame>(
              Instant const&)> intrinsic_acceleration,
      RelativeDegreesOfFreedom const& tolerance,
      AdaptiveStepSizeIntegrator<TimedBurnMotion> integrator,
      Instant const& t);

 private:
  std::vector<not_null<std::unique_ptr<MassiveBody>>> bodies_;
  std::map<not_null<MassiveBody const*>,
           ContinuousTrajectory<Frame>> trajectories_;

  FixedStepSizeIntegrator<PlanetaryMotion> const& planetary_integrator_;
};

}  // namespace physics
}  // namespace principia
