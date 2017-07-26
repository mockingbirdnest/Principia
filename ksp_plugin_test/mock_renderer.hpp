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
  MOCK_METHOD1(SetPlottingFrameConstRef,
               void(NavigationFrame const& plotting_frame));

  MOCK_CONST_METHOD0(GetPlottingFrame, not_null<NavigationFrame const*> ());

  not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderBarycentricTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation)
      const override;
  MOCK_CONST_METHOD6(
      FillRenderedBarycentricTrajectoryInWorld,
      void(Instant const& time,
           DiscreteTrajectory<Barycentric>::Iterator const& begin,
           DiscreteTrajectory<Barycentric>::Iterator const& end,
           Position<World> const& sun_world_position,
           Rotation<Barycentric, AliceSun> const& planetarium_rotation,
           std::unique_ptr<DiscreteTrajectory<World>>*
               rendered_barycentric_trajectory_in_world));

  MOCK_CONST_METHOD1(
      BarycentricToWorldSun,
      OrthogonalMap<Barycentric, WorldSun>(
          Rotation<Barycentric, AliceSun> const& planetarium_rotation));
  MOCK_CONST_METHOD3(
      PlottingToWorld,
      RigidTransformation<Navigation, World>(
          Instant const& time,
          Position<World> const& sun_world_position,
          Rotation<Barycentric, AliceSun> const& planetarium_rotation));
  MOCK_CONST_METHOD3(
      WorldToPlotting,
      RigidTransformation<World, Navigation>(
          Instant const& time,
          Position<World> const& sun_world_position,
          Rotation<Barycentric, AliceSun> const& planetarium_rotation));
};

}  // namespace internal_renderer

using internal_renderer::MockRenderer;

}  // namespace ksp_plugin
}  // namespace principia
