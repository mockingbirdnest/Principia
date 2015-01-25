#pragma once

#include "physics/massless_body.hpp"

namespace principia {
namespace physics {

inline bool MasslessBody::is_massless() const {
  return true;
}

inline bool MasslessBody::is_oblate() const {
  return false;
}

inline not_null<std::unique_ptr<MasslessBody>> MasslessBody::ReadFromMessage(
    serialization::Body const& message) {
  CHECK(message.HasExtension(serialization::MasslessBody::massless_body));
  return ReadFromMessage(
      message.GetExtension(serialization::MasslessBody::massless_body));
}

not_null<std::unique_ptr<MasslessBody>> MasslessBody::ReadFromMessage(
      serialization::MasslessBody const& message) {
  return std::make_unique<MasslessBody>();
}

}  // namespace physics
}  // namespace principia
