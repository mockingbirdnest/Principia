#pragma once

#include "physics/body.hpp"

#include "physics/oblate_body.hpp"

namespace principia {
namespace physics {

template<typename Frame>
bool Body::is_compatible_with() const {
  return CompatibilityHelper<Frame, 
                             Frame::is_inertial>::is_compatible_with(this);
}

template<typename Frame>
class Body::CompatibilityHelper<Frame, false> {
   public:
     static bool is_compatible_with(Body const* body);
};

template<typename Frame>
class Body::CompatibilityHelper<Frame, true> {
   public:
     static bool is_compatible_with(Body const* body);
};

template<typename Frame>
bool Body::CompatibilityHelper<Frame, false>::is_compatible_with(
    Body const* body) {
  return !body->is_oblate();
}

template<typename Frame>
bool Body::CompatibilityHelper<Frame, true>::is_compatible_with(
    Body const* body) {
  return !body->is_oblate() ||
         dynamic_cast<OblateBody<Frame> const*>(body) != nullptr;
}

}  // namespace physics
}  // namespace principia
