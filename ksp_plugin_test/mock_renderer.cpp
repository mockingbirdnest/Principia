#include "mock_renderer.hpp"

#include "physics/mock_dynamic_frame.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using physics::MockDynamicFrame;

MockCelestial* const sun = new MockCelestial;

MockRenderer::MockRenderer()
    : Renderer(sun,
               std::make_unique<MockDynamicFrame<Barycentric, Navigation>>()) {}

void MockRenderer::SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) {
  SetPlottingFrameConstRef(*plotting_frame);
  plotting_frame.release();
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
MockRenderer::RenderBarycentricTrajectoryInWorld(
    Instant const& time,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  std::unique_ptr<DiscreteTrajectory<World>>
      rendered_barycentric_trajectory_in_world;
  FillRenderedBarycentricTrajectoryInWorld(
      time,
      begin,
      end,
      sun_world_position,
      planetarium_rotation,
      &rendered_barycentric_trajectory_in_world);
  return std::move(rendered_barycentric_trajectory_in_world);
}

}  // namespace internal_renderer
}  // namespace ksp_plugin
}  // namespace principia
