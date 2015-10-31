#include "ksp_plugin/mock_plugin.hpp"

#include <vector>

#include "gmock/gmock.h"

namespace principia {
namespace ksp_plugin {

MockPlugin::MockPlugin()
    : Plugin(Instant(),
      Angle()) {}

void MockPlugin::DirectlyInsertCelestial(
    Index const celestial_index,
    Index const* const parent_index,
    DegreesOfFreedom<Barycentric> const & initial_state,
    std::unique_ptr<MassiveBody> body) {
  DirectlyInsertCelestialConstRef(celestial_index,
                                  parent_index,
                                  initial_state,
                                  body);
}

not_null<std::unique_ptr<RenderingFrame>>
MockPlugin::NewBodyCentredNonRotatingRenderingFrame(
    Index const reference_body_index) const {
  std::unique_ptr<RenderingFrame> rendering_frame;
  FillBodyCentredNonRotatingRenderingFrame(reference_body_index,
                                           &rendering_frame);
  return std::move(rendering_frame);
}

not_null<std::unique_ptr<RenderingFrame>>
MockPlugin::NewBarycentricRotatingRenderingFrame(
    Index const primary_index,
    Index const secondary_index) const {
  std::unique_ptr<RenderingFrame> rendering_frame;
  FillBarycentricRotatingRenderingFrame(primary_index,
                                        secondary_index,
                                        &rendering_frame);
  return std::move(rendering_frame);
}

void MockPlugin::AddVesselToNextPhysicsBubble(
    GUID const& vessel_guid,
    std::vector<IdAndOwnedPart> parts) {
  AddVesselToNextPhysicsBubbleConstRef(vessel_guid, parts);
}

}  // namespace ksp_plugin
}  // namespace principia
