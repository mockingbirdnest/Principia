#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"

namespace principia {

namespace physics {

template<typename Frame>
FrameField<Frame> CoordinateFrame() {
  auto const identically_identity =
      [](Position<Frame> const&) -> Rotation<Frame, Frame> {
        return Rotation<Frame, Frame>::Identity();
      };
  return identically_identity;
}

}  // namespace physics
}  // namespace principia
