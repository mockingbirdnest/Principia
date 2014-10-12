#pragma once

#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
Celestial<Frame>::Celestial(
    GravitationalParameter const& gravitational_parameter)
    : body(new Body<Frame>(gravitational_parameter)) {}

}  // namespace ksp_plugin
}  // namespace principia
