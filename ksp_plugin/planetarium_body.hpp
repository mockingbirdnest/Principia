#pragma once

#include "ksp_plugin/planetarium.hpp"

#include <algorithm>

#include "base/algebra.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/sign.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/hermite3.hpp"
#include "numerics/quadrature.hpp"
#include "physics/similar_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace _planetarium {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_sign;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_quadrature;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

// A helper function that converts a trajectory expressed in `Frame` into one
// expressed in `Navigation`, using the given `plotting_frame` if needed for the
// transformation.
template<typename Frame>
DegreesOfFreedom<Navigation> EvaluateDegreesOfFreedomInNavigation(
    PlottingFrame const& plotting_frame,
    Trajectory<Frame> const& trajectory,
    Instant const& t);

template<>
inline DegreesOfFreedom<Navigation>
EvaluateDegreesOfFreedomInNavigation<Barycentric>(
    PlottingFrame const& plotting_frame,
    Trajectory<Barycentric> const& trajectory,
    Instant const& t) {
  SimilarMotion<Barycentric, Navigation> to_plotting_frame_at_t =
      plotting_frame.ToThisFrameAtTimeSimilarly(t);
  return to_plotting_frame_at_t(trajectory.EvaluateDegreesOfFreedom(t));
}

template<>
inline DegreesOfFreedom<Navigation>
EvaluateDegreesOfFreedomInNavigation<Navigation>(
    PlottingFrame const& plotting_frame,
    Trajectory<Navigation> const& trajectory,
    Instant const& t) {
  return trajectory.EvaluateDegreesOfFreedom(t);
}

template<typename Frame>
void Planetarium::PlotMethod3(
    Trajectory<Frame> const& trajectory,
    Instant const& first_time,
    Instant const& last_time,
    bool const reverse,
    std::function<void(ScaledSpacePoint const&)> const& add_point,
    int const max_points,
    Length* const minimal_distance) const {
  double const tan²_angular_resolution =
      Pow<2>(parameters_.tan_angular_resolution_);
  auto const final_time = reverse ? first_time : last_time;
  auto previous_time = reverse ? last_time : first_time;

  if (minimal_distance != nullptr) {
    *minimal_distance = Infinity<Length>;
  }

  Sign const direction = reverse ? Sign::Negative() : Sign::Positive();
  if (direction * (final_time - previous_time) <= Time{}) {
    return;
  }
  DegreesOfFreedom<Navigation> const initial_degrees_of_freedom =
      EvaluateDegreesOfFreedomInNavigation<Frame>(
          *plotting_frame_, trajectory, previous_time);
  Position<Navigation> previous_position =
      initial_degrees_of_freedom.position();
  Velocity<Navigation> previous_velocity =
      initial_degrees_of_freedom.velocity();
  Time Δt = final_time - previous_time;

  add_point(plotting_to_scaled_space_(previous_time, previous_position));
  int points_added = 1;

  Instant t;
  double estimated_tan²_error;
  std::optional<DegreesOfFreedom<Navigation>> degrees_of_freedom;
  Position<Navigation> position;
  Square<Length> minimal_squared_distance = Infinity<Square<Length>>;

  goto estimate_tan²_error;

  while (points_added < max_points &&
         direction * (previous_time - final_time) < Time{}) {
    do {
      // One square root because we have squared errors, another one because the
      // errors are quadratic in time (in other words, two square roots because
      // the squared errors are quartic in time).
      // A safety factor prevents catastrophic retries.
      Δt *= 0.9 * Sqrt(Sqrt(tan²_angular_resolution / estimated_tan²_error));
    estimate_tan²_error:
      t = previous_time + Δt;
      if (direction * (t - final_time) > Time{}) {
        t = final_time;
        Δt = t - previous_time;
      }
      Position<Navigation> const extrapolated_position =
          previous_position + previous_velocity * Δt;
      degrees_of_freedom = EvaluateDegreesOfFreedomInNavigation<Frame>(
          *plotting_frame_, trajectory, t);
      position = degrees_of_freedom->position();

      // The quadratic term of the error between the linear interpolation and
      // the actual function is maximized halfway through the segment, so it is
      // 1/2 (Δt/2)² f″(t-Δt) = (1/2 Δt² f″(t-Δt)) / 4; the squared error is
      // thus (1/2 Δt² f″(t-Δt))² / 16.
      estimated_tan²_error =
          perspective_.Tan²AngularDistance(extrapolated_position, position) /
          16;
    } while (estimated_tan²_error > tan²_angular_resolution);

    previous_time = t;
    previous_position = position;
    previous_velocity = degrees_of_freedom->velocity();

    add_point(plotting_to_scaled_space_(t, position));
    ++points_added;

    if (minimal_distance != nullptr) {
      minimal_squared_distance =
          std::min(minimal_squared_distance,
                   perspective_.SquaredDistanceFromCamera(position));
    }
  }
  if (minimal_distance != nullptr) {
    *minimal_distance = Sqrt(minimal_squared_distance);
  }
}

template<typename Frame>
void Planetarium::PlotMethod4(
    Trajectory<Frame> const& trajectory,
    Instant const& first_time,
    Instant const& last_time,
    bool const reverse,
    std::function<void(ScaledSpacePoint const&)> const& add_point,
    int const max_points,
    Length* const minimal_distance) const {
  auto const final_time = reverse ? first_time : last_time;
  auto previous_time = reverse ? last_time : first_time;

  if (minimal_distance != nullptr) {
    *minimal_distance = Infinity<Length>;
  }

  Sign const direction = reverse ? Sign::Negative() : Sign::Positive();
  if (direction * (final_time - previous_time) <= Time{}) {
    return;
  }
  DegreesOfFreedom<Navigation> const initial_degrees_of_freedom =
      EvaluateDegreesOfFreedomInNavigation<Frame>(
          *plotting_frame_, trajectory, previous_time);
  Position<Navigation> previous_position =
      initial_degrees_of_freedom.position();
  Vector<AngularFrequency, Navigation> previous_projected_velocity =
      ProperMotion(initial_degrees_of_freedom) *
          Normalize(previous_position - perspective_.camera());

  Time Δt = final_time - previous_time;

  add_point(plotting_to_scaled_space_(previous_time, previous_position));
  int points_added = 1;

  Instant t;
  Angle rms_apparent_distance;
  Position<Navigation> position;
  Vector<AngularFrequency, Navigation> projected_velocity;
  Square<Length> minimal_squared_distance = Infinity<Square<Length>>;

  goto estimate_tan²_error;

  while (points_added < max_points &&
         direction * (previous_time - final_time) < Time{}) {
    do {
      // Consider the case where the trajectory is an arc of circle of angle 2 α
      // and the segment is the corresponding chord.  The higher-order term of
      // the integral below is 4 α⁵ / 15, which yields an RMS of O(α²) after
      // division by the time interval and application of the square root.
      // Therefore, our errors are of degree 2 in time.
      // A safety factor prevents catastrophic retries.
      Δt *= 0.9 * Sqrt(parameters_.angular_resolution_ / rms_apparent_distance);
    estimate_tan²_error:
      t = previous_time + Δt;
      if (direction * (t - final_time) > Time{}) {
        t = final_time;
        Δt = t - previous_time;
      }

      auto const degrees_of_freedom =
          EvaluateDegreesOfFreedomInNavigation<Frame>(
              *plotting_frame_, trajectory, t);
      position = degrees_of_freedom.position();

      // The velocity on the linear segment from the previous position to the
      // current one.
      Velocity<Navigation> const linear_velocity =
          (position - previous_position) / Δt;

      // We now compute a 3rd-degree Hermite approximation of the displacement
      // between the trajectory and the line segment linking the previous and
      // current positions. The displacement is zero at the extremities.  The
      // projected velocities are obtained by applying the rotation vector to
      // the actual velocities and are represented as a (somewhat unusual)
      // `Vector<AngularFrequency, Navigation>`.  Note that neither the line
      // segment nor the trajectory are on the sphere, but for sufficiently
      // small time intervals they are close
      projected_velocity =
          ProperMotion(degrees_of_freedom) *
          Normalize(position - perspective_.camera());
      auto const previous_projected_linear_velocity =
          ProperMotion({previous_position, linear_velocity}) *
          Normalize(previous_position - perspective_.camera());
      auto const projected_linear_velocity =
          ProperMotion({position, linear_velocity}) *
          Normalize(position - perspective_.camera());

      Hermite3<Vector<Angle, Navigation>, Instant> const error_approximation(
          {previous_time, t},
          {Vector<Angle, Navigation>{}, Vector<Angle, Navigation>{}},
          {previous_projected_velocity - previous_projected_linear_velocity,
           projected_velocity - projected_linear_velocity});

      // Our metric is the root mean square of the norm of the displacement
      // between the trajectory and the linear segment, normalized by the
      // duration of the time interval.  It is an angle.  Note that the square
      // norm is a polynomial of degree 6, and the 4-point Gauss-Legendre method
      // is exact for polynomials of degree 7.
      rms_apparent_distance =
          Sqrt(GaussLegendre<4>(
                   [&error_approximation](Instant const& time) {
                     return error_approximation.Evaluate(time).Norm²();
                   },
                   previous_time,
                   t) /
               Δt);
    } while (rms_apparent_distance > parameters_.angular_resolution_);

    previous_time = t;
    previous_position = position;
    previous_projected_velocity = projected_velocity;

    add_point(plotting_to_scaled_space_(t, position));
    ++points_added;

    if (minimal_distance != nullptr) {
      minimal_squared_distance =
          std::min(minimal_squared_distance,
                   perspective_.SquaredDistanceFromCamera(position));
    }
  }
  if (minimal_distance != nullptr) {
    *minimal_distance = Sqrt(minimal_squared_distance);
  }
}

inline AngularVelocity<Navigation> Planetarium::ProperMotion(
    DegreesOfFreedom<Navigation> const& degrees_of_freedom) const {
  auto const r = degrees_of_freedom.position() - perspective_.camera();
  return Wedge(r, degrees_of_freedom.velocity()) / r.Norm²() * Radian;
}

}  // namespace internal
}  // namespace _planetarium
}  // namespace ksp_plugin
}  // namespace principia
