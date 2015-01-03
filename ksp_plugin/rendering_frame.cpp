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
  Transforms<Barycentric, Rendering, Barycentric>::
      LazyTrajectory<Barycentric> const body_prolongation =
          std::bind(&Celestial::prolongation, &body_);
  auto transforms(
      Transforms<Barycentric, Rendering, Barycentric>::BodyCentredNonRotating(
          body_prolongation,
          body_prolongation));
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
  Transforms<Barycentric, Rendering, Barycentric>::
      LazyTrajectory<Barycentric> const primary_prolongation =
          std::bind(&Celestial::prolongation, &primary_);
  Transforms<Barycentric, Rendering, Barycentric>::
      LazyTrajectory<Barycentric> const secondary_prolongation =
          std::bind(&Celestial::prolongation, &secondary_);
  auto transforms(
      Transforms<Barycentric, Rendering, Barycentric>::BarycentricRotating(
          primary_prolongation,
          primary_prolongation,
          secondary_prolongation,
          secondary_prolongation));
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
