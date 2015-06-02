#pragma once

#include "physics/ephemeris.hpp"

#include <functional>

#include "base/map_util.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "physics/continuous_trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using base::FindOrDie;
using geometry::InnerProduct;
using geometry::R3Element;
using integrators::IntegrationProblem;
using quantities::Exponentiation;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace physics {

namespace {

template<typename B>
std::enable_if_t<std::is_base_of<Body, B>::value, not_null<B const*>>
ConvertTo(not_null<Body const*> const body) {
// Dynamic casting is expensive, as in 3x slower for the benchmarks.  Do that in
// debug mode to catch bugs, but not in optimized mode where we want all the
// performance we can get.
#ifdef _DEBUG
  return dynamic_cast<B const*>(static_cast<Body const*>(body));
#else
  return static_cast<not_null<B const*>>(body);
#endif
}

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

    auto const inserted = bodies_to_trajectories_.emplace(
                              std::piecewise_construct,
                              std::forward_as_tuple(body.get()),
                              std::forward_as_tuple(step_,
                                                    low_fitting_tolerance_,
                                                    high_fitting_tolerance_));
    CHECK(inserted.second);
    ContinuousTrajectory<Frame>* const trajectory = &inserted.first->second;

    if (body->is_oblate()) {
      // Inserting at the beginning of the vectors is O(N).
      oblate_bodies_.insert(oblate_bodies_.begin(), body.get());
      bodies_.insert(bodies_.begin(), std::move(body));
      oblate_trajectories_.insert(oblate_trajectories_.begin(), trajectory);
      last_state_.positions.insert(last_state_.positions.begin(),
                                   degrees_of_freedom.position());
      last_state_.velocities.insert(last_state_.velocities.begin(),
                                    degrees_of_freedom.velocity());
      ++number_of_oblate_bodies_;
    } else {
      // Inserting at the end of the vectors is O(1).
      spherical_bodies_.push_back(body.get());
      bodies_.push_back(std::move(body));
      spherical_trajectories_.push_back(trajectory);
      last_state_.positions.push_back(degrees_of_freedom.position());
      last_state_.velocities.push_back(degrees_of_freedom.velocity());
      ++number_of_spherical_bodies_;
    }
  }

  equation_.compute_acceleration =
      std::bind(&Ephemeris::ComputeGravitationalAccelerations,
                this, _1, _2, _3);
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
  problem.append_state = std::bind(&Ephemeris::AppendState, this, _1);
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
  for (auto& trajectory : oblate_trajectories_) {
    trajectory->Append(
        state.time.value,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value));
    ++index;
  }
  for (auto& trajectory : spherical_trajectories_) {
    trajectory->Append(
        state.time.value,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value));
    ++index;
  }
}

template<typename Frame>
template<bool body1_is_oblate,
         bool body2_is_oblate>
inline void Ephemeris<Frame>::ComputeOneBodyGravitationalAcceleration(
    MassiveBody const& body1,
    size_t const b1,
    std::vector<not_null<MassiveBody const*>> const& bodies2,
    size_t const b2_begin,
    size_t const b2_end,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations) {
  GravitationalParameter const& μ1 =
      body1.gravitational_parameter();
  for (std::size_t b2 = std::max(b1 + 1, b2_begin); b2 < b2_end; ++b2) {
    MassiveBody const& body2 = *bodies2[b2 - b2_begin];
    GravitationalParameter const& μ2 =
        body2.gravitational_parameter();

    Displacement<Frame> const Δq = positions[b1] - positions[b2];

    Exponentiation<Length, 2> const Δq_squared = InnerProduct(Δq, Δq);
    // NOTE(phl): Don't try to compute one_over_Δq_squared here, it makes the
    // non-oblate path slower.
    Exponentiation<Length, -3> const one_over_Δq_cubed =
        Sqrt(Δq_squared) / (Δq_squared * Δq_squared);

    auto const μ1_over_Δq_cubed = μ1 * one_over_Δq_cubed;
    (*accelerations)[b2] += Δq * μ1_over_Δq_cubed;

    // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
    // sive corporum duorum actiones in se mutuo semper esse æquales &
    // in partes contrarias dirigi.
    auto const μ2_over_Δq_cubed = μ2 * one_over_Δq_cubed;
    (*accelerations)[b1] -= Δq * μ2_over_Δq_cubed;

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
      }
      if (body2_is_oblate) {
        Vector<Acceleration, Frame> const order_2_zonal_acceleration2 =
            Order2ZonalAcceleration<Frame>(
                static_cast<OblateBody<Frame> const&>(body2),
                Δq,
                one_over_Δq_squared,
                one_over_Δq_cubed);
        (*accelerations)[b1] -= order_2_zonal_acceleration2;
      }
    }
  }
}

template<typename Frame>
void Ephemeris<Frame>::ComputeGravitationalAccelerations(
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations) {
  accelerations->assign(accelerations->size(), Vector<Acceleration, Frame>());

  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *oblate_bodies_[b1];
    ComputeOneBodyGravitationalAcceleration<true /*body1_is_oblate*/,
                                            true /*body2_is_oblate*/>(
        body1, b1,
        oblate_bodies_ /*bodies2*/,
        0 /*b2_begin*/,
        number_of_oblate_bodies_ /*b2_end*/,
        positions,
        accelerations);
    ComputeOneBodyGravitationalAcceleration<true /*body1_is_oblate*/,
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
    ComputeOneBodyGravitationalAcceleration<false /*body1_is_oblate*/,
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

}  // namespace physics
}  // namespace principia
