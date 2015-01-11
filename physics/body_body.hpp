#pragma once

#include "physics/body.hpp"

#include "physics/oblate_body.hpp"

using principia::base::check_not_null;

namespace principia {
namespace physics {

template<typename Frame>
bool Body::is_compatible_with() const {
  return CompatibilityHelper<Frame,
                             Frame::is_inertial>::is_compatible_with(check_not_null(this));
}

template<typename Frame>
class Body::CompatibilityHelper<Frame, false> {
 public:
  static bool is_compatible_with(not_null<Body const*> const body);
};

template<typename Frame>
class Body::CompatibilityHelper<Frame, true> {
 public:
  static bool is_compatible_with(not_null<Body const*> const body);
};

template<typename Frame>
bool Body::CompatibilityHelper<Frame, false>::is_compatible_with(
    not_null<Body const*> const body) {
  return !body->is_oblate();
}

template<typename Frame>
bool Body::CompatibilityHelper<Frame, true>::is_compatible_with(
    not_null<Body const*> const body) {
  return !body->is_oblate() ||
         dynamic_cast<OblateBody<Frame> const*>(
            static_cast<Body const*>(body)) != nullptr;
}

}  // namespace physics
}  // namespace principia
