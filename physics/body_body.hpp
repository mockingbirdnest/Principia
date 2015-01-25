#pragma once

#include "physics/body.hpp"

#include "physics/oblate_body.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"

namespace principia {
namespace physics {

template<typename Frame>
bool Body::is_compatible_with() const {
  return CompatibilityHelper<Frame,
                             Frame::is_inertial>::is_compatible_with(this);
}

inline not_null<std::unique_ptr<Body>> Body::ReadFromMessage(
    serialization::Body const& message) {
  if (message.HasExtension(serialization::MasslessBody::massless_body)) {
    return MasslessBody::ReadFromMessage(
        message.GetExtension(serialization::MasslessBody::massless_body));
  } else if (message.HasExtension(serialization::MassiveBody::massive_body)) {
    return MassiveBody::ReadFromMessage(
        message.GetExtension(serialization::MassiveBody::massive_body));
  } else {
    LOG(FATAL) << "serialization::Body is neither massive nor massless";
    base::noreturn();
  }
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
