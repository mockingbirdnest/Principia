#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/trajectory.hpp"

namespace principia {

using geometry::Position;
using geometry::Vector;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::FixedStepSizeIntegrator;
using integrators::SpecialSecondOrderDifferentialEquation;

namespace physics {

template<typename Frame>
class Ephemeris {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  // The equation describing the motion of the |bodies_|.
  using NewtonianMotionEquation =
      SpecialSecondOrderDifferentialEquation<Position<Frame>>;

  // Constructs an Ephemeris that owns the |bodies|.  The elements of vectors
  // |bodies| and |initial_state| correspond to one another.
  Ephemeris(std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies,
            std::vector<DegreesOfFreedom<Frame>> const& initial_state,
            Instant const& initial_time,
            FixedStepSizeIntegrator<NewtonianMotionEquation> const&
                planetary_integrator,
            Time const& step,
            Length const& low_fitting_tolerance,
            Length const& high_fitting_tolerance);

  ContinuousTrajectory<Frame> const& trajectory(
      not_null<MassiveBody const*> body) const;

  // The maximum of the |t_min|s of the trajectories.
  Instant t_min() const;
  // The mimimum of the |t_max|s of the trajectories.
  Instant t_max() const;

  // Calls |ForgetBefore| on all trajectories.
  void ForgetBefore(Instant const& t);

  // Prolongs the ephemeris up to at least |t|.  After the call, |t_max() >= t|.
  void Prolong(Instant const& t);

  // Integrates, until exactly |t|, the |trajectory| followed by a massless body
  // in the gravitational potential described by |*this|, and subject to the
  // given |intrinsic_acceleration|.
  // If |t > t_max()|, calls |Prolong(t)| beforehand.
  // The |length_| and |speed_integration_tolerance|s are used to compute the
  // |tolerance_to_error_ratio| for step size control.
  // TODO(phl): Remove intrinsic_acceleration?  It is in the trajectory.
  void Flow(not_null<Trajectory<Frame>*> const trajectory,
            Length const& length_integration_tolerance,
            Speed const& speed_integration_tolerance,
            AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const&
                integrator,
            Instant const& t);

 private:
  void AppendMassiveBodiesState(
      typename NewtonianMotionEquation::SystemState const& state);
  static void AppendMasslessBodyState(
      typename NewtonianMotionEquation::SystemState const& state,
      not_null<Trajectory<Frame>*> const trajectory);

  // Computes the acceleration due to one body, |body1| (with index |b1| in the
  // |positions| and |accelerations| arrays) on the bodies |bodies2| (with
  // indices [b2_begin, b2_end[ in the |positions| and |accelerations| arrays).
  // The template parameters specify what we know about the bodies, and
  // therefore what forces apply.
  template<bool body1_is_oblate,
           bool body2_is_oblate>
  static void ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies(
      MassiveBody const& body1,
      size_t const b1,
      std::vector<not_null<MassiveBody const*>> const& bodies2,
      size_t const b2_begin,
      size_t const b2_end,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations);

  // Computes the acceleration due to one body, |body1| (with index |b1| in the
  // |hints|, |bodies_| and |trajectories_| arrays) on a massless body at the
  // given |position|.  The template parameter specifies what we know about the
  // massive body, and therefore what forces apply.
  template<bool body1_is_oblate>
  void ComputeGravitationalAccelerationByMassiveBodyOnMasslessBody(
      Instant const& t,
      MassiveBody const& body1,
      size_t const b1,
      Position<Frame> const& position,
      not_null<Vector<Acceleration, Frame>*> const acceleration,
      not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
          const hints);

  // Computes the accelerations between all the massive bodies in |bodies_|.
  void ComputeMassiveBodiesGravitationalAccelerations(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations);

  // Computes the acceleration exerted by the massive bodies in |bodies_| on a
  // massless body.  The massless body may have an intrinsic acceleration
  // described in its |trajectory| object.  The |hints| are passed to
  // ComputeGravitationalAccelerationByMassiveBodyOnMasslessBody for efficient
  // computation of the positions of the massive bodies.
  void ComputeMasslessBodyGravitationalAccelerations(
      not_null<Trajectory<Frame> const*> const trajectory,
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
      not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
          const hints);

  // Computes an estimate of the ratio |tolerance / error|.
  static double ToleranceToErrorRatio(
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance,
      Time const& current_step_size,
      typename NewtonianMotionEquation::SystemStateError const& error);

  // The oblate bodies precede the spherical bodies in this vector.  The system
  // state is indexed in the same order.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies_;

  // The indices in |bodies_| correspond to those in |oblate_bodies_| and
  // |spherical_bodies_|, in sequence.  The elements of |oblate_bodies_| are
  // really |OblateBody<Frame>| but it's inconvenient to express.
  std::vector<not_null<MassiveBody const*>> oblate_bodies_;
  std::vector<not_null<MassiveBody const*>> spherical_bodies_;

  // The indices in |bodies_| correspond to those in |trajectories_|.
  std::vector<not_null<ContinuousTrajectory<Frame>*>> trajectories_;

  std::map<not_null<MassiveBody const*>, ContinuousTrajectory<Frame>>
      bodies_to_trajectories_;

  // This will refer to a static object returned by a factory.
  FixedStepSizeIntegrator<NewtonianMotionEquation> const& planetary_integrator_;
  Time const step_;
  Length const low_fitting_tolerance_;
  Length const high_fitting_tolerance_;
  typename NewtonianMotionEquation::SystemState last_state_;

  int number_of_spherical_bodies_ = 0;
  int number_of_oblate_bodies_ = 0;

  NewtonianMotionEquation massive_bodies_equation_;
};

}  // namespace physics
}  // namespace principia

#include "physics/ephemeris_body.hpp"
