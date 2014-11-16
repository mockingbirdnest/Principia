#include "ksp_plugin/rendering_frame.hpp"

#include <utility>

#include "physics/transforms.hpp"

using principia::physics::Transforms;

namespace principia {
namespace ksp_plugin {

BodyCentredNonRotatingFrame::BodyCentredNonRotatingFrame(
    Celestial const& body) : body_(body) {}

std::unique_ptr<Trajectory<Barycentric>>
BodyCentredNonRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentric> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentric>> result =
      std::make_unique<Trajectory<Barycentric>>(actual_trajectory.body());
  // TODO(phl): Should tag the frames differently.
  auto transforms(
      Transforms<Barycentric, Barycentric, Barycentric>::BodyCentredNonRotating(
          body_.prolongation(),
          body_.prolongation()));
  auto actual_it = transforms.first(&actual_trajectory);
  auto body_it = body_.prolongation().on_or_after(actual_it.time());

  // First build the trajectory resulting from the first transform.
  Trajectory<Barycentric> intermediate_trajectory(actual_trajectory.body());
  for (; !actual_it.at_end(); ++actual_it, ++body_it) {
    // Advance over the bits of the actual trajectory that don't have a matching
    // time in the body trajectory.
    while (actual_it.time() != body_it.time()) {
      ++actual_it;
    }
    intermediate_trajectory.Append(actual_it.time(),
                                   actual_it.degrees_of_freedom());
  }

  // Then build the final trajectory using the second transform.
  for (auto intermediate_it = transforms.second(&intermediate_trajectory);
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
      std::make_unique<Trajectory<Barycentric>>(actual_trajectory.body());
  // TODO(phl): Should tag the frames differently.
  auto transforms(
      Transforms<Barycentric, Barycentric, Barycentric>::BarycentricRotating(
          primary_.prolongation(),
          primary_.prolongation(),
          secondary_.prolongation(),
          secondary_.prolongation()));
  auto actual_it = transforms.first(&actual_trajectory);
  auto primary_it = primary_.prolongation().on_or_after(actual_it.time());
  auto secondary_it = secondary_.prolongation().on_or_after(actual_it.time());

  // First build the trajectory resulting from the first transform.
  Trajectory<Barycentric> intermediate_trajectory(actual_trajectory.body());
  for (; !actual_it.at_end(); ++actual_it, ++primary_it, ++secondary_it) {
    // Advance over the bits of the actual trajectory that don't have a matching
    // time in the trajectories.
    while (actual_it.time() != primary_it.time() ||
           actual_it.time() != secondary_it.time()) {
      ++actual_it;
    }
    intermediate_trajectory.Append(actual_it.time(),
                                   actual_it.degrees_of_freedom());
  }

  // Then build the final trajectory using the second transform.
  for (auto intermediate_it = transforms.second(&intermediate_trajectory);
       !intermediate_it.at_end();
       ++intermediate_it) {
    result->Append(intermediate_it.time(),
                   intermediate_it.degrees_of_freedom());
  }

  return std::move(result);
}

}  // namespace ksp_plugin
}  // namespace principia
