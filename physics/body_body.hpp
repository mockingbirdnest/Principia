
#pragma once

#include "physics/body.hpp"

#include "base/not_constructible.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"

namespace principia {
namespace physics {
namespace internal_body {

using base::dynamic_cast_not_null;
using base::not_constructible;

template<typename Frame>
bool Body::is_compatible_with() const {
  return CompatibilityHelper<Frame,
                             Frame::is_inertial>::is_compatible_with(this);
}

inline not_null<std::unique_ptr<Body>> Body::ReadFromMessage(
    serialization::Body const& message) {
  if (message.has_massless_body()) {
    return MasslessBody::ReadFromMessage(message.massless_body());
  } else if (message.has_massive_body()) {
    return MassiveBody::ReadFromMessage(message.massive_body());
  } else {
    LOG(FATAL) << "Body is neither massive nor massless";
    base::noreturn();
  }
}

template<typename Frame>
struct Body::CompatibilityHelper<Frame, false> : not_constructible {
  static bool is_compatible_with(not_null<Body const*> body);
};

template<typename Frame>
struct Body::CompatibilityHelper<Frame, true> : not_constructible {
  static bool is_compatible_with(not_null<Body const*> body);
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
         dynamic_cast_not_null<OblateBody<Frame> const*>(body) != nullptr;
}

}  // namespace internal_body
}  // namespace physics
}  // namespace principia
