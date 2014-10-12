#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
Vessel<Frame>::Vessel(Celestial<Frame> const* parent)
    : body(new Body<Frame>(GravitationalParameter())),
      parent(CHECK_NOTNULL(parent)) {}

}  // namespace ksp_plugin
}  // namespace principia
