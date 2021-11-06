#pragma once

#include "ksp_plugin/renderer.hpp"

#include "gmock/gmock.h"
#include "physics/mock_dynamic_frame.hpp"
#include "ksp_plugin_test/mock_celestial.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using physics::MockDynamicFrame;

MockCelestial* const sun = new MockCelestial;

class MockRenderer : public Renderer {
 public:
  MockRenderer()
      : Renderer(
            sun,
            std::make_unique<MockDynamicFrame<Barycentric, Navigation>>()){};

  MOCK_METHOD(void,
              SetPlottingFrame,
              (not_null<std::unique_ptr<NavigationFrame>> plotting_frame),
              (override));

  MOCK_METHOD(not_null<NavigationFrame const*>,
              GetPlottingFrame,
              (),
              (const, override));

  MOCK_METHOD(DiscreteTraject0ry<World>,
              RenderBarycentricTrajectoryInWorld,
              (Instant const& time,
               DiscreteTraject0ry<Barycentric>::iterator const& begin,
               DiscreteTraject0ry<Barycentric>::iterator const& end,
               Position<World> const& sun_world_position,
               (Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
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
