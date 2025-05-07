#pragma once

#include "ksp_plugin/planetarium.hpp"

#include <algorithm>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/pair.hpp"
#include "geometry/sign.hpp"
#include "numerics/hermite3.hpp"
#include "numerics/quadrature.hpp"
#include "physics/similar_motion.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _planetarium {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_pair;
using namespace principia::geometry::_sign;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_quadrature;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

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

  Vector<Inverse<Time>, Navigation> previous_projected_velocity;
  std::tie(std::ignore, previous_projected_velocity) =
      SphericalProjection(initial_degrees_of_freedom);

  Time Δt = final_time - previous_time;

  add_point(plotting_to_scaled_space_(previous_time, previous_position));
  int points_added = 1;

  Instant t;
  double metric;
  Position<Navigation> position;
  Vector<Inverse<Time>, Navigation> projected_velocity;
  Velocity<Navigation> velocity;
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
      Δt *= 0.9 * Sqrt(tan²_angular_resolution / metric);
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
      velocity = degrees_of_freedom.velocity();

      Position<Navigation> const linear_midpoint =
          Barycentre({previous_position, position});
      Velocity<Navigation> const linear_velocity =
          (position - previous_position) / Δt;

      std::tie(std::ignore, projected_velocity) =
          SphericalProjection(degrees_of_freedom);
      auto const& [_1, previous_projected_linear_velocity] =
          SphericalProjection({previous_position, linear_velocity});
      auto const& [_2, projected_linear_velocity] =
          SphericalProjection({position, linear_velocity});

      Hermite3<Vector<double, Navigation>, Instant> const error_approximation(
          {previous_time, t},
          {Vector<double, Navigation>{}, Vector<double, Navigation>{}},
          {previous_projected_velocity - previous_projected_linear_velocity,
           projected_velocity - projected_linear_velocity});

      metric = Sqrt(GaussLegendre<4>(
                        [&error_approximation](Instant const& time) {
                          return error_approximation.Evaluate(time).Norm²();
                        },
                        previous_time,
                        t) /
                    Δt);
    } while (metric > tan²_angular_resolution);

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

inline std::pair<Point<Vector<double, Navigation>>,
                 Vector<Inverse<Time>, Navigation>>
Planetarium::SphericalProjection(
    DegreesOfFreedom<Navigation> const& degrees_of_freedom) const {
  auto const r = degrees_of_freedom.position() - perspective_.camera();
  return {Point<Vector<double, Navigation>>{} + Normalize(r),
          degrees_of_freedom.velocity().OrthogonalizationAgainst(r) / r.Norm()};
}

}  // namespace internal
}  // namespace _planetarium
}  // namespace ksp_plugin
}  // namespace principia
