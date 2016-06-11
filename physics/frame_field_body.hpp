
#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_frame_field {

template<typename Frame>
FrameField<Frame> CoordinateFrame() {
  return [](Position<Frame> const&) -> Rotation<Frame, Frame> {
           return Rotation<Frame, Frame>::Identity();
         };
}

}  // namespace internal_frame_field
}  // namespace physics
}  // namespace principia
