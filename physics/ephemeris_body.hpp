#pragma once

#include "physics/ephemeris.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "base/map_util.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "physics/continuous_trajectory.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using base::FindOrDie;
using geometry::InnerProduct;
using geometry::R3Element;
using integrators::AdaptiveStepSize;
using integrators::IntegrationProblem;
using quantities::Abs;
using quantities::Exponentiation;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace physics {

namespace {

// If j is a unit vector along the axis of rotation, and r is the separation
// between the bodies, the acceleration computed here is:
//
//   -(J2 / |r|^5) (3 j (r.j) + r (3 - 15 (r.j)^2 / |r|^2) / 2)
//
// Where |r| is the norm of r and r.j is the inner product.
template<typename Frame>
FORCE_INLINE Vector<Acceleration, Frame>
    Order2ZonalAcceleration(
        OblateBody<Frame> const& body,
        Displacement<Frame> const& r,
        Exponentiation<Length, -2> const& one_over_r_squared,
        Exponentiation<Length, -3> const& one_over_r_cubed) {
  Vector<double, Frame> const& axis = body.axis();
  Length const r_axis_projection = InnerProduct(axis, r);
  auto const j2_over_r_fifth =
      body.j2() * one_over_r_cubed * one_over_r_squared;
  Vector<Acceleration, Frame> const& axis_acceleration =
      (-3 * j2_over_r_fifth * r_axis_projection) * axis;
  Vector<Acceleration, Frame> const& radial_acceleration =
      (j2_over_r_fifth *
           (-1.5 +
            7.5 * r_axis_projection *
                  r_axis_projection * one_over_r_squared)) * r;
  return axis_acceleration + radial_acceleration;
}

}  // namespace

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies,
    std::vector<DegreesOfFreedom<Frame>> const& initial_state,
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

  for (int i = 0; i < bodies.size(); ++i) {
    auto& body = bodies[i];
    DegreesOfFreedom<Frame> const& degrees_of_freedom = initial_state[i];

    unowned_bodies_.push_back(body.get());

    auto const inserted = bodies_to_trajectories_.emplace(
                              std::piecewise_construct,
                              std::forward_as_tuple(body.get()),
                              std::forward_as_tuple(step_,
                                                    low_fitting_tolerance_,
                                                    high_fitting_tolerance_));
    CHECK(inserted.second);
    ContinuousTrajectory<Frame>* const trajectory = &inserted.first->second;
    trajectory->Append(initial_time, degrees_of_freedom);

    if (body->is_oblate()) {
      // Inserting at the beginning of the vectors is O(N).
      oblate_bodies_.insert(oblate_bodies_.begin(), body.get());
      bodies_.insert(bodies_.begin(), std::move(body));
      trajectories_.insert(trajectories_.begin(), trajectory);
      last_state_.positions.insert(last_state_.positions.begin(),
                                   degrees_of_freedom.position());
      last_state_.velocities.insert(last_state_.velocities.begin(),
                                    degrees_of_freedom.velocity());
      ++number_of_oblate_bodies_;
    } else {
      // Inserting at the end of the vectors is O(1).
      spherical_bodies_.push_back(body.get());
      bodies_.push_back(std::move(body));
      trajectories_.push_back(trajectory);
      last_state_.positions.push_back(degrees_of_freedom.position());
      last_state_.velocities.push_back(degrees_of_freedom.velocity());
      ++number_of_spherical_bodies_;
    }
  }

  massive_bodies_equation_.compute_acceleration =
      std::bind(&Ephemeris::ComputeMassiveBodiesGravitationalAccelerations,
                this, _1, _2, _3);
}

template<typename Frame>
std::vector<MassiveBody const*> const& Ephemeris<Frame>::bodies() const {
  return unowned_bodies_;
}

template<typename Frame>
ContinuousTrajectory<Frame> const& Ephemeris<Frame>::trajectory(
    not_null<MassiveBody const*> body) const {
  return FindOrDie(bodies_to_trajectories_, body);
}

template<typename Frame>
bool Ephemeris<Frame>::empty() const {
  for (auto const& pair : bodies_to_trajectories_) {
    ContinuousTrajectory<Frame> const& trajectory = pair.second;
    if (trajectory.empty()) {
      return true;
    }
  }
  return false;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_min() const {
  Instant t_min;
  for (auto const& pair : bodies_to_trajectories_) {
    ContinuousTrajectory<Frame> const& trajectory = pair.second;
    t_min = std::max(t_min, trajectory.t_min());
  }
  return t_min;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_max() const {
  Instant t_max = bodies_to_trajectories_.begin()->second.t_max();
  for (auto const& pair : bodies_to_trajectories_) {
    ContinuousTrajectory<Frame> const& trajectory = pair.second;
    t_max = std::min(t_max, trajectory.t_max());
  }
  return t_max;
}

template<typename Frame>
void Ephemeris<Frame>::ForgetBefore(Instant const& t) {
  for (auto const& pair : bodies_to_trajectories_) {
    ContinuousTrajectory<Frame>& trajectory = pair.second;
    trajectory.ForgetBefore(t);
  }
}

template<typename Frame>
void Ephemeris<Frame>::Prolong(Instant const& t) {
  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation = massive_bodies_equation_;
  problem.append_state =
      std::bind(&Ephemeris::AppendMassiveBodiesState, this, _1);
  problem.t_final = t;
  problem.initial_state = &last_state_;

  // Perform the integration.  Note that we may have to iterate until |t_max()|
  // actually reaches |t| because the last series may not be fully determined
  // after the first integration.
  do {
    planetary_integrator_.Solve(problem, step_);
    // Here |problem.initial_state| still points at |last_state_|, which is the
    // state at the end of the previous call to |Solve|.  It is therefore the
    // right initial state for the next call to |Solve|, if any.
    problem.t_final += step_;
  } while (t_max() < t);
}

template<typename Frame>
void Ephemeris<Frame>::FlowWithAdaptiveStep(
    not_null<Trajectory<Frame>*> const trajectory,
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance,
    AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
    Instant const& t) {
  if (empty() || t > t_max()) {
    Prolong(t);
  }

  std::vector<typename ContinuousTrajectory<Frame>::Hint> hints(bodies_.size());
  NewtonianMotionEquation massless_body_equation;
  massless_body_equation.compute_acceleration =
      std::bind(&Ephemeris::ComputeMasslessBodyGravitationalAccelerations,
                this, trajectory, _1, _2, _3, &hints);

  typename NewtonianMotionEquation::SystemState initial_state;
  auto const trajectory_last = trajectory->last();
  auto const last_degrees_of_freedom = trajectory_last.degrees_of_freedom();
  initial_state.time = trajectory_last.time();
  initial_state.positions.push_back(last_degrees_of_freedom.position());
  initial_state.velocities.push_back(last_degrees_of_freedom.velocity());

  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation = massless_body_equation;
  problem.append_state =
      std::bind(&Ephemeris::AppendMasslessBodyState, _1, trajectory);
  problem.t_final = t;
  problem.initial_state = &initial_state;

  AdaptiveStepSize<NewtonianMotionEquation> step_size;
  step_size.first_time_step = problem.t_final - initial_state.time.value;
  step_size.safety_factor = 0.9;
  step_size.tolerance_to_error_ratio =
      std::bind(&Ephemeris<Frame>::ToleranceToErrorRatio,
                std::cref(length_integration_tolerance),
                std::cref(speed_integration_tolerance),
                _1, _2);

  integrator.Solve(problem, step_size);
}

template<typename Frame>
void Ephemeris<Frame>::FlowWithFixedStep(
    std::vector<not_null<Trajectory<Frame>*>> const& trajectories,
    Time const& step,
    Instant const& t) {
  if (empty() || t > t_max()) {
    Prolong(t);
  }

  std::vector<typename ContinuousTrajectory<Frame>::Hint> hints(bodies_.size());
  NewtonianMotionEquation massless_body_equation;
  massless_body_equation.compute_acceleration =
      std::bind(&Ephemeris::ComputeMasslessBodyGravitationalAccelerations,
                this, trajectory, _1, _2, _3, &hints);

  typename NewtonianMotionEquation::SystemState initial_state;
    for (auto const& trajectory : trajectories) {
    auto const trajectory_last = trajectory->last();
    auto const last_degrees_of_freedom = trajectory_last.degrees_of_freedom();
    initial_state.time = trajectory_last.time();
    initial_state.positions.push_back(last_degrees_of_freedom.position());
    initial_state.velocities.push_back(last_degrees_of_freedom.velocity());
  }

  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation = massless_body_equation;
  problem.append_state =
      std::bind(&Ephemeris::AppendMasslessBodiesState, _1, trajectories;
  problem.t_final = t;
  problem.initial_state = &initial_state;

  planetary_integrator_.Solve(problem, step);
}

template<typename Frame>
void Ephemeris<Frame>::AppendMassiveBodiesState(
    typename NewtonianMotionEquation::SystemState const& state) {
  last_state_ = state;
  int index = 0;
  for (auto& trajectory : trajectories_) {
    trajectory->Append(
        state.time.value,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value));
    ++index;
  }
}

template<typename Frame>
void Ephemeris<Frame>::AppendMasslessBodyState(
      typename NewtonianMotionEquation::SystemState const& state,
      not_null<Trajectory<Frame>*> const trajectory) {
  trajectory->Append(state.time.value,
                     DegreesOfFreedom<Frame>(state.positions[0].value,
                                             state.velocities[0].value));
}

template<typename Frame>
template<bool body1_is_oblate,
         bool body2_is_oblate>
void Ephemeris<Frame>::
ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies(
    MassiveBody const& body1,
    size_t const b1,
    std::vector<not_null<MassiveBody const*>> const& bodies2,
    size_t const b2_begin,
    size_t const b2_end,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations) {
  Vector<Acceleration, Frame>& acceleration_on_b1 = (*accelerations)[b1];
  GravitationalParameter const& μ1 = body1.gravitational_parameter();
  for (std::size_t b2 = std::max(b1 + 1, b2_begin); b2 < b2_end; ++b2) {
    Vector<Acceleration, Frame>& acceleration_on_b2 = (*accelerations)[b2];
    MassiveBody const& body2 = *bodies2[b2 - b2_begin];
    GravitationalParameter const& μ2 = body2.gravitational_parameter();

    Displacement<Frame> const Δq = positions[b1] - positions[b2];

    Exponentiation<Length, 2> const Δq_squared = InnerProduct(Δq, Δq);
    // NOTE(phl): Don't try to compute one_over_Δq_squared here, it makes the
    // non-oblate path slower.
    Exponentiation<Length, -3> const one_over_Δq_cubed =
        Sqrt(Δq_squared) / (Δq_squared * Δq_squared);

    auto const μ1_over_Δq_cubed = μ1 * one_over_Δq_cubed;
    acceleration_on_b2 += Δq * μ1_over_Δq_cubed;

    // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
    // sive corporum duorum actiones in se mutuo semper esse æquales &
    // in partes contrarias dirigi.
    auto const μ2_over_Δq_cubed = μ2 * one_over_Δq_cubed;
    acceleration_on_b1 -= Δq * μ2_over_Δq_cubed;

    if (body1_is_oblate || body2_is_oblate) {
      Exponentiation<Length, -2> const one_over_Δq_squared = 1 / Δq_squared;
      if (body1_is_oblate) {
        Vector<Acceleration, Frame> const order_2_zonal_acceleration1 =
            Order2ZonalAcceleration<Frame>(
                static_cast<OblateBody<Frame> const &>(body1),
                Δq,
                one_over_Δq_squared,
                one_over_Δq_cubed);
        (*accelerations)[b2] += order_2_zonal_acceleration1;
        // Don't update |(*accelerations)[b1]| here: we exempt ourselves from
        // Newton's third law because the oblate bodies are generally much
        // bigger than the other bodies, so the ratio of the masses makes the
        // opposite acceleration very small.
      }
      if (body2_is_oblate) {
        Vector<Acceleration, Frame> const order_2_zonal_acceleration2 =
            Order2ZonalAcceleration<Frame>(
                static_cast<OblateBody<Frame> const&>(body2),
                Δq,
                one_over_Δq_squared,
                one_over_Δq_cubed);
        (*accelerations)[b1] -= order_2_zonal_acceleration2;
        // Yet another exemption from Newton's third law.
      }
    }
  }
}

template<typename Frame>
template<bool body1_is_oblate>
void Ephemeris<Frame>::
ComputeGravitationalAccelerationByMassiveBodyOnMasslessBody(
    Instant const& t,
    MassiveBody const& body1,
    size_t const b1,
    Position<Frame> const& position,
    not_null<Vector<Acceleration, Frame>*> const acceleration,
    not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
        const hints) {
  GravitationalParameter const& μ1 = body1.gravitational_parameter();

  Displacement<Frame> const Δq =
      trajectories_[b1]->EvaluatePosition(t, &(*hints)[b1]) - position;

  Exponentiation<Length, 2> const Δq_squared = InnerProduct(Δq, Δq);
  // NOTE(phl): Don't try to compute one_over_Δq_squared here, it makes the
  // non-oblate path slower.
  Exponentiation<Length, -3> const one_over_Δq_cubed =
      Sqrt(Δq_squared) / (Δq_squared * Δq_squared);

  auto const μ1_over_Δq_cubed = μ1 * one_over_Δq_cubed;
  *acceleration += Δq * μ1_over_Δq_cubed;

  if (body1_is_oblate) {
    Exponentiation<Length, -2> const one_over_Δq_squared = 1 / Δq_squared;
    Vector<Acceleration, Frame> const order_2_zonal_acceleration1 =
        Order2ZonalAcceleration<Frame>(
            static_cast<OblateBody<Frame> const &>(body1),
            Δq,
            one_over_Δq_squared,
            one_over_Δq_cubed);
    *acceleration += order_2_zonal_acceleration1;
  }
}

template<typename Frame>
void Ephemeris<Frame>::ComputeMassiveBodiesGravitationalAccelerations(
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations) {
  accelerations->assign(accelerations->size(), Vector<Acceleration, Frame>());

  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *oblate_bodies_[b1];
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        true /*body1_is_oblate*/,
        true /*body2_is_oblate*/>(
        body1, b1,
        oblate_bodies_ /*bodies2*/,
        0 /*b2_begin*/,
        number_of_oblate_bodies_ /*b2_end*/,
        positions,
        accelerations);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        true /*body1_is_oblate*/,
        false /*body2_is_oblate*/>(
        body1, b1,
        spherical_bodies_ /*bodies2*/,
        number_of_oblate_bodies_ /*b2_begin*/,
        number_of_oblate_bodies_ +
            number_of_spherical_bodies_ /*b2_end*/,
        positions,
        accelerations);
  }
  for (std::size_t b1 = number_of_oblate_bodies_;
       b1 < number_of_oblate_bodies_ +
            number_of_spherical_bodies_;
       ++b1) {
    MassiveBody const& body1 =
        *spherical_bodies_[b1 - number_of_oblate_bodies_];
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        false /*body1_is_oblate*/,
        false /*body2_is_oblate*/>(
        body1, b1,
        spherical_bodies_ /*bodies2*/,
        number_of_oblate_bodies_ /*b2_begin*/,
        number_of_oblate_bodies_ +
            number_of_spherical_bodies_ /*b2_end*/,
        positions,
        accelerations);
  }
}

template<typename Frame>
void Ephemeris<Frame>::ComputeMasslessBodyGravitationalAccelerations(
      not_null<Trajectory<Frame> const*> const trajectory,
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
      not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
          const hints) {
  CHECK_EQ(1, positions.size());
  CHECK_EQ(1, accelerations->size());
  Position<Frame> const& position = positions[0];
  Vector<Acceleration, Frame>& acceleration = (*accelerations)[0];
  acceleration = Vector<Acceleration, Frame>();

  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *oblate_bodies_[b1];
    ComputeGravitationalAccelerationByMassiveBodyOnMasslessBody<
        true /*body1_is_oblate*/>(
        t,
        body1, b1,
        position,
        &acceleration,
        hints);
  }
  for (std::size_t b1 = number_of_oblate_bodies_;
       b1 < number_of_oblate_bodies_ +
            number_of_spherical_bodies_;
       ++b1) {
    MassiveBody const& body1 =
        *spherical_bodies_[b1 - number_of_oblate_bodies_];
    ComputeGravitationalAccelerationByMassiveBodyOnMasslessBody<
        false /*body1_is_oblate*/>(
        t,
        body1, b1,
        position,
        &acceleration,
        hints);
  }
  // Finally, take into account the intrinsic accelerations.
  if (trajectory->has_intrinsic_acceleration()) {
    Vector<Acceleration, Frame> const intrinsic_acceleration =
        trajectory->evaluate_intrinsic_acceleration(t);
    acceleration += intrinsic_acceleration;
  }
}

template<typename Frame>
double Ephemeris<Frame>::ToleranceToErrorRatio(
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance,
    Time const& current_step_size,
    typename NewtonianMotionEquation::SystemStateError const& error) {
  Length max_length_error;
  Speed max_speed_error;
  for (auto const& position_error : error.position_error) {
    max_length_error = std::max(max_length_error,
                                position_error.Norm());
  }
  for (auto const& velocity_error : error.velocity_error) {
    max_speed_error = std::max(max_speed_error,
                               velocity_error.Norm());
  }
  return std::min(length_integration_tolerance / max_length_error,
                  speed_integration_tolerance / max_speed_error);
}

}  // namespace physics
}  // namespace principia
