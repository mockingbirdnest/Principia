#include "ksp_plugin/mock_plugin.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace ksp_plugin {

MockPlugin::MockPlugin()
    : Plugin(Instant(),
      Index(0),
      1 * SIUnit<GravitationalParameter>(),
      Angle()) {}

std::unique_ptr<BodyCentredNonRotatingFrame>
MockPlugin::NewBodyCentredNonRotatingFrame(
    Index const reference_body_index) const {
  std::unique_ptr<BodyCentredNonRotatingFrame> frame;
  FillBodyCentredNonRotatingFrame(reference_body_index, &frame);
  return frame;
}

std::unique_ptr<BarycentricRotatingFrame>
MockPlugin::NewBarycentricRotatingFrame(
    Index const primary_index,
    Index const secondary_index) const {
  std::unique_ptr<BarycentricRotatingFrame> frame;
  FillBarycentricRotatingFrame(primary_index, secondary_index, &frame);
  return frame;
}

}  // namespace ksp_plugin
}  // namespace principia
