
#include "ksp_plugin_test/mock_plugin.hpp"

#include <vector>

#include "gmock/gmock.h"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

MockPlugin::MockPlugin() : Plugin(Instant(), Instant(), Angle()) {}

void MockPlugin::InsertCelestialAbsoluteCartesian(
      Index const celestial_index,
      std::experimental::optional<Index> const& parent_index,
      DegreesOfFreedom<Barycentric> const& initial_state,
      base::not_null<std::unique_ptr<RotatingBody<Barycentric> const>> body) {
  InsertCelestialAbsoluteCartesianConstRef(
      celestial_index, parent_index, initial_state, body);
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
MockPlugin::RenderBarycentricTrajectoryInWorld(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position) const {
  std::unique_ptr<DiscreteTrajectory<World>>
      rendered_barycentric_trajectory_in_world;
  FillRenderedBarycentricTrajectoryInWorld(
      begin,
      end,
      sun_world_position,
      &rendered_barycentric_trajectory_in_world);
  return std::move(rendered_barycentric_trajectory_in_world);
}

not_null<std::unique_ptr<NavigationFrame>>
MockPlugin::NewBodyCentredNonRotatingNavigationFrame(
    Index const reference_body_index) const {
  std::unique_ptr<NavigationFrame> navigation_frame;
  FillBodyCentredNonRotatingNavigationFrame(reference_body_index,
                                            &navigation_frame);
  return std::move(navigation_frame);
}

not_null<std::unique_ptr<NavigationFrame>>
MockPlugin::NewBarycentricRotatingNavigationFrame(
    Index const primary_index,
    Index const secondary_index) const {
  std::unique_ptr<NavigationFrame> navigation_frame;
  FillBarycentricRotatingNavigationFrame(primary_index,
                                         secondary_index,
                                         &navigation_frame);
  return std::move(navigation_frame);
}

void MockPlugin::SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) {
  SetPlottingFrameConstRef(*plotting_frame);
  plotting_frame.release();
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
