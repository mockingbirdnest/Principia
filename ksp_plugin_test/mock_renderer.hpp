#pragma once

#include <memory>

#include "gmock/gmock.h"
#include "ksp_plugin_test/mock_celestial.hpp"  // 🧙 For MockCelestial.
#include "physics/mock_rigid_reference_frame.hpp"  // 🧙 For MockRigidReferenceFrame.  // NOLINT
#include "physics/rigid_reference_frame.hpp"

namespace principia {
namespace ksp_plugin {
namespace _renderer {
namespace internal {

using namespace principia::physics::_rigid_reference_frame;

MockCelestial* const sun = new MockCelestial;

class MockRenderer : public Renderer {
 public:
  MockRenderer()
      : Renderer(sun,
                 std::make_unique<
                     MockRigidReferenceFrame<Barycentric, Navigation>>()){};

  MOCK_METHOD(void,
              SetPlottingFrame,
              (not_null<std::unique_ptr<PlottingFrame>> plotting_frame),
              (override));

  MOCK_METHOD(not_null<PlottingFrame const*>,
              GetPlottingFrame,
              (),
              (const, override));

  MOCK_METHOD(DiscreteTrajectory<World>,
              RenderBarycentricTrajectoryInWorld,
              (Instant const& time,
               DiscreteTrajectory<Barycentric>::iterator const& begin,
               DiscreteTrajectory<Barycentric>::iterator const& end,
               Position<World> const& sun_world_position,
               (Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
              (const, override));

  MOCK_METHOD((OrthogonalMap<Barycentric, WorldSun>),
              BarycentricToWorldSun,
              ((Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
              (const, override));
  MOCK_METHOD((Similarity<Navigation, World>),
              PlottingToWorld,
              (Instant const& time,
               Position<World> const& sun_world_position,
               (Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
              (const, override));
  MOCK_METHOD((Similarity<World, Navigation>),
              WorldToPlotting,
              (Instant const& time,
               Position<World> const& sun_world_position,
               (Rotation<Barycentric, AliceSun> const& planetarium_rotation)),
              (const, override));
};

}  // namespace internal

using internal::MockRenderer;

}  // namespace _renderer
}  // namespace ksp_plugin
}  // namespace principia
