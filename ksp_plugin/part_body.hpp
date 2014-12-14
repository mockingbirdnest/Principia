#pragma once

#include "part.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
inline Part<Frame>::Part(
    DegreesOfFreedom<Frame> const& degrees_of_freedom,
    Mass const& mass,
    Vector<Acceleration, Frame> const&
        gravitational_acceleration_to_be_applied_by_ksp)
    : degrees_of_freedom(degrees_of_freedom),
      mass(mass),
      gravitational_acceleration_to_be_applied_by_ksp(
          gravitational_acceleration_to_be_applied_by_ksp) {}

template<typename Frame>
std::ostream& operator<<(std::ostream& out, Part<Frame> const& part) {
  return out << "{"
      << part.degrees_of_freedom << ", "
      << part.mass << ", "
      << part.gravitational_acceleration_to_be_applied_by_ksp << "}";
}

}  // namespace ksp_plugin
}  // namespace principia
