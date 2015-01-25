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

inline void MasslessBody::WriteToMessage(
    not_null<serialization::Body*> message) const {
  WriteToMessage(message->mutable_massless_body());
}

inline void MasslessBody::WriteToMessage(
    not_null<serialization::MasslessBody*> message) const {}

inline not_null<std::unique_ptr<MasslessBody>> MasslessBody::ReadFromMessage(
    serialization::Body const& message) {
  CHECK(message.has_massless_body());
  return ReadFromMessage(message.massless_body());
}

inline not_null<std::unique_ptr<MasslessBody>> MasslessBody::ReadFromMessage(
    serialization::MasslessBody const& message) {
  return std::make_unique<MasslessBody>();
}

}  // namespace physics
}  // namespace principia
