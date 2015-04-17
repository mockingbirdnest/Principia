#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"

namespace principia {

namespace physics {

template<typename Frame>
FrameField<Frame> CoordinateField() {
  auto const identity =
      [](Position<Frame>, Instant) {
        return Rotation<Frame, Frame>::Identity();
      };
  return identity;
}

}  // namespace physics
}  // namespace principia
