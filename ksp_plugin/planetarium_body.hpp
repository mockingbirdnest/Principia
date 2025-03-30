#pragma once

#include "ksp_plugin/planetarium.hpp"

#include <algorithm>

#include "geometry/sign.hpp"
#include "physics/similar_motion.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _planetarium {
namespace internal {

using namespace principia::geometry::_sign;
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
    double const tan_angular_resolution,
    bool const reverse,
    std::function<void(ScaledSpacePoint const&)> const& add_point,
    int const max_points,
    Length* const minimal_distance) const {
  //LOG(INFO)<<"rev="<<reverse<<" first="<<first_time<<" last="<<last_time;
  auto const final_time = reverse ? first_time : last_time;
  auto previous_time = reverse ? last_time : first_time;

  if (minimal_distance != nullptr) {
    *minimal_distance = Infinity<Length>;
  }

  Sign const direction = reverse ? Sign::Negative() : Sign::Positive();
  if (direction * (final_time - previous_time) <= Time{}) {
    return;
  }
///Position only
  DegreesOfFreedom<Navigation> const initial_degrees_of_freedom =
      EvaluateDegreesOfFreedomInNavigation<Frame>(
          *plotting_frame_, trajectory, previous_time);
  Position<Navigation> previous_position =
      initial_degrees_of_freedom.position();
  Time Δt = (final_time - previous_time) / 3;

  add_point(plotting_to_scaled_space_(previous_position));
  int points_added = 1;

  Instant t;
  double estimated_tan_error;
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
///Comment
      // A safety factor prevents catastrophic retries.
      Δt *= 0.9 * Sqrt(tan_angular_resolution / estimated_tan_error);
    estimate_tan²_error:
      t = previous_time + Δt;
      if (direction * (t - final_time) > Time{}) {
        t = final_time;
        Δt = (t - previous_time) / 2;
      }
      degrees_of_freedom = EvaluateDegreesOfFreedomInNavigation<Frame>(
          *plotting_frame_, trajectory, t);
      position = degrees_of_freedom->position();

      Instant next_t = t + Δt;
      if (direction * (next_t - final_time) > Time{}) {
        next_t = final_time;
      }
      DegreesOfFreedom<Navigation> const next_degrees_of_freedom =
          EvaluateDegreesOfFreedomInNavigation<Frame>(
              *plotting_frame_, trajectory, next_t);
      Position<Navigation> const next_position =
          next_degrees_of_freedom.position();

      ///TanAngleBetween?
      Angle const α = AngleBetween(next_position - position,
                                   position - previous_position);
      estimated_tan_error = Abs(Tan(α));

//      LOG_IF(INFO, !reverse)<<"previous="<<previous_position;
//      LOG_IF(INFO, !reverse)<<"current="<<position;
//      LOG_IF(INFO, !reverse)<<"next="<<next_position;
//      LOG_IF(INFO, !reverse)<<"t= "<<t<<" Δt= "<<Δt<<
//" error="<<estimated_tan_error<<" α="<<α;

    } while (estimated_tan_error > tan_angular_resolution);

    previous_time = t;
    previous_position = position;

    add_point(plotting_to_scaled_space_(position));
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

}  // namespace internal
}  // namespace _planetarium
}  // namespace ksp_plugin
}  // namespace principia
