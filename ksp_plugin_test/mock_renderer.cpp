#include "mock_renderer.hpp"

#include "physics/mock_dynamic_frame.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using physics::MassiveBody;
using physics::MockDynamicFrame;
using physics::RotatingBody;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

RotatingBody<Barycentric> body(MassiveBody::Parameters(1 * Kilogram),
                               RotatingBody<Barycentric>::Parameters(
                                   /*mean_radius=*/1 * Metre,
                                   /*reference_angle=*/0 * Radian,
                                   /*reference_instant=*/Instant{},
                                   /*angular_frequency=*/10 * Radian / Second,
                                   /*ascension_of_pole=*/0 * Radian,
                                   /*declination_of_pole=*/π / 2 * Radian));
Celestial const sun(&body);

MockRenderer::MockRenderer()
    : Renderer(&sun,
               std::make_unique<MockDynamicFrame<Barycentric, Navigation>>()) {}

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
