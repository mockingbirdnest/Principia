#pragma once

#include "part.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
inline Part<Frame>::Part(
    DegreesOfFreedom<Frame> const& degrees_of_freedom,
    Mass const& mass,
    Vector<Acceleration, Frame> const& expected_ksp_gravity)
    : degrees_of_freedom(degrees_of_freedom),
      mass(mass),
      expected_ksp_gravity(expected_ksp_gravity) {}

}  // namespace ksp_plugin
}  // namespace principia
