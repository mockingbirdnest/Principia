#pragma once

#include "ksp_plugin/renderer.hpp"

#include "gmock/gmock.h"
#include "ksp_plugin_test/mock_celestial.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

class MockRenderer : public Renderer {
 public:
  MockRenderer();

  // NOTE(phl): Needed because gMock wants to copy the unique_ptr<>.
  void SetPlottingFrame(
      not_null<std::unique_ptr<NavigationFrame>> plotting_frame) override;
  MOCK_METHOD(void,
              SetPlottingFrameConstRef,
              (NavigationFrame const& plotting_frame),
              (override));

  MOCK_METHOD(not_null<NavigationFrame const*>,
              GetPlottingFrame,
              (),
              (const, override));

  not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderBarycentricTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation)
      const override;
  MOCK_METHOD(void,
              FillRenderedBarycentricTrajectoryInWorld,
              (Instant const& time,
               DiscreteTrajectory<Barycentric>::Iterator const& begin,
               DiscreteTrajectory<Barycentric>::Iterator const& end,
               Position<World> const& sun_world_position,
               (Rotation<Barycentric, AliceSun> const& planetarium_rotation),
               std::unique_ptr<DiscreteTrajectory<World>>*
                   rendered_barycentric_trajectory_in_world),
              (const, override));

  MOCK_METHOD((OrthogonalMap<Barycentric, WorldSun>),
              BarycentricToWorldSun,
              ((Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
              (const, override));
  MOCK_METHOD((RigidTransformation<Navigation, World>),
              PlottingToWorld,
              (Instant const& time,
               Position<World> const& sun_world_position,
               (Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
              (const, override));
  MOCK_METHOD((RigidTransformation<World, Navigation>),
              WorldToPlotting,
              (Instant const& time,
               Position<World> const& sun_world_position,
               (Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
              (const, override));
};

}  // namespace internal_renderer

using internal_renderer::MockRenderer;

}  // namespace ksp_plugin
}  // namespace principia
