#pragma once

#include "ksp_plugin/planetarium.hpp"

#include <algorithm>

#include "geometry/sign.hpp"
#include "numerics/hermite3.hpp"
#include "numerics/quadrature.hpp"
#include "physics/similar_motion.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace _planetarium {
namespace internal {

using namespace principia::geometry::_sign;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_quadrature;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_elementary_functions;
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
      // One square root because we have squared errors, another one because the
      // errors are quadratic in time (in other words, two square roots because
      // the squared errors are quartic in time).
      // A safety factor prevents catastrophic retries.
      //TODO(phl)Comment
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
      // between the trajectory and the linear segment linking the previous and
      // current positions. The displacement is zero at the extremities.  The
      // projected velocities are obtained by applying the rotation vector to
      // the actual velocities and are represented as a (somewhat unusual)
      // `Vector<AngularFrequency, Navigation>`.  Note that neither the linear
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
      // duration of the time interval.  It is an angle.
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
