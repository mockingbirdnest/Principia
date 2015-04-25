#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"

namespace principia {

namespace physics {

template<typename Frame>
FrameField<Frame> CoordinateFrame() {
  return [](Position<Frame> const&) -> Rotation<Frame, Frame> {
           return Rotation<Frame, Frame>::Identity();
         };
}

}  // namespace physics
}  // namespace principia
