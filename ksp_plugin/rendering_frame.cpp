#include "ksp_plugin/rendering_frame.hpp"

#include <utility>

#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "physics/body.hpp"
#include "physics/transforms.hpp"

using principia::geometry::Vector;
using principia::physics::Body;
using principia::physics::Transforms;

namespace principia {
namespace ksp_plugin {

BodyCentredNonRotatingFrame::BodyCentredNonRotatingFrame(
    Celestial const& body) : body_(body) {}

std::unique_ptr<Trajectory<Barycentric>>
BodyCentredNonRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentric> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentric>> result =
      std::make_unique<Trajectory<Barycentric>>(actual_trajectory.body<Body>());
  auto transforms(
      Transforms<Barycentric, Rendering, Barycentric>::BodyCentredNonRotating(
          body_.prolongation(),
          body_.prolongation()));
  auto actual_it = transforms->first(&actual_trajectory);

  // First build the trajectory resulting from the first transform.
  Trajectory<Rendering> intermediate_trajectory(actual_trajectory.body<Body>());
  for (; !actual_it.at_end(); ++actual_it) {
    intermediate_trajectory.Append(actual_it.time(),
                                   actual_it.degrees_of_freedom());
  }

  // Then build the final trajectory using the second transform.
  for (auto intermediate_it = transforms->second(&intermediate_trajectory);
       !intermediate_it.at_end();
       ++intermediate_it) {
    result->Append(intermediate_it.time(),
                   intermediate_it.degrees_of_freedom());
  }

  return std::move(result);
}

BarycentricRotatingFrame::BarycentricRotatingFrame(
    Celestial const& primary,
    Celestial const& secondary)
    : primary_(primary),
      secondary_(secondary) {}

std::unique_ptr<Trajectory<Barycentric>>
BarycentricRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentric> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentric>> result =
      std::make_unique<Trajectory<Barycentric>>(actual_trajectory.body<Body>());
  auto transforms(
      Transforms<Barycentric, Rendering, Barycentric>::BarycentricRotating(
          primary_.prolongation(),
          primary_.prolongation(),
          secondary_.prolongation(),
          secondary_.prolongation()));
  auto actual_it = transforms->first(&actual_trajectory);

  // First build the trajectory resulting from the first transform.
  Trajectory<Rendering> intermediate_trajectory(actual_trajectory.body<Body>());
  for (; !actual_it.at_end(); ++actual_it) {
    intermediate_trajectory.Append(actual_it.time(),
                                   actual_it.degrees_of_freedom());
  }

  // Then build the final trajectory using the second transform.
  for (auto intermediate_it = transforms->second(&intermediate_trajectory);
       !intermediate_it.at_end();
       ++intermediate_it) {
    result->Append(intermediate_it.time(),
                   intermediate_it.degrees_of_freedom());
  }

  return std::move(result);
}

}  // namespace ksp_plugin
}  // namespace principia
