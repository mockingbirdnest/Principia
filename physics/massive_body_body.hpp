#pragma once

#include "physics/massive_body.hpp"

#include "geometry/frame.hpp"
#include "glog/logging.h"
#include "physics/oblate_body.hpp"
#include "quantities/constants.hpp"

using principia::constants::GravitationalConstant;
using principia::geometry::UnknownInertialFrame;

namespace principia {
namespace physics {

inline MassiveBody::MassiveBody(
    GravitationalParameter const& gravitational_parameter)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant) {
  CHECK_NE(gravitational_parameter, GravitationalParameter())
      << "Massive body cannot have zero gravitational parameter";
}

inline MassiveBody::MassiveBody(Mass const& mass)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass) {
  CHECK_NE(mass, Mass())
      << "Massive body cannot have zero mass";
}

inline GravitationalParameter const&
MassiveBody::gravitational_parameter() const {
  return gravitational_parameter_;
}

inline Mass const& MassiveBody::mass() const {
  return mass_;
}

inline bool MassiveBody::is_massless() const {
  return false;
}

inline bool MassiveBody::is_oblate() const {
  return false;
}

inline void MassiveBody::WriteToMessage(
    not_null<serialization::Body*> const message) const {
  WriteToMessage(message->mutable_massive_body());
}

inline void MassiveBody::WriteToMessage(
    not_null<serialization::MassiveBody*> const message) const {
  gravitational_parameter_.WriteToMessage(
      message->mutable_gravitational_parameter());
}

inline not_null<std::unique_ptr<MassiveBody>> MassiveBody::ReadFromMessage(
    serialization::Body const& message) {
  CHECK(message.has_massive_body());
  return ReadFromMessage(message.massive_body());
}

inline not_null<std::unique_ptr<MassiveBody>> MassiveBody::ReadFromMessage(
    serialization::MassiveBody const& message) {
  if (message.HasExtension(serialization::OblateBody::oblate_body)) {
    return OblateBody<UnknownInertialFrame>::ReadFromMessage(message);
  } else {
    return std::make_unique<MassiveBody>(
        GravitationalParameter::ReadFromMessage(
            message.gravitational_parameter()));
  }
}

}  // namespace physics
}  // namespace principia
