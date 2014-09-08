#pragma once

#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {

template<typename Frame>
DegreesOfFreedom<Frame>::DegreesOfFreedom(
    Point<Vector<Length, Frame>> const& position,
    Vector<Speed, Frame> const& velocity)
    : position(position),
      velocity(velocity) {}

}  // namespace physics
}  // namespace principia
