#pragma once

#include "part.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
inline Part<Frame>::Part(
    Position<Frame> const& position,
    Velocity<Frame> const& velocity,
    Mass const& mass,
    Vector<Acceleration, Frame> const& expected_ksp_gravity)
    : position(position),
      velocity(velocity),
      mass(mass),
      expected_ksp_gravity(expected_ksp_gravity) {}


}  // namespace ksp_plugin
}  // namespace principia
