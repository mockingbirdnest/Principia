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

}  // namespace internal_renderer
}  // namespace ksp_plugin
}  // namespace principia
