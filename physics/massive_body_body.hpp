#pragma once

#include "physics/massive_body.hpp"

#include <string>

#include "base/macros.hpp"
#include "geometry/frame.hpp"
#include "glog/logging.h"
#include "google/protobuf/descriptor.h"
#include "physics/oblate_body.hpp"
#include "quantities/constants.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using quantities::constants::GravitationalConstant;
using geometry::Frame;
using geometry::ReadFrameFromMessage;

namespace physics {

inline MassiveBody::Parameters::Parameters(
    GravitationalParameter const& gravitational_parameter)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant) {
  CHECK_NE(gravitational_parameter, GravitationalParameter())
      << "Massive body cannot have zero gravitational parameter";
}

inline MassiveBody::Parameters::Parameters(Mass const& mass)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass) {
  CHECK_NE(mass, Mass()) << "Massive body cannot have zero mass";
}

inline MassiveBody::MassiveBody(Parameters const& parameters)
    : parameters_(parameters) {}

inline GravitationalParameter const&
MassiveBody::gravitational_parameter() const {
  return parameters_.gravitational_parameter_;
}

inline Mass const& MassiveBody::mass() const {
  return parameters_.mass_;
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
  parameters_.gravitational_parameter_.WriteToMessage(
      message->mutable_gravitational_parameter());
}

inline not_null<std::unique_ptr<MassiveBody>> MassiveBody::ReadFromMessage(
    serialization::Body const& message) {
  CHECK(message.has_massive_body());
  return ReadFromMessage(message.massive_body());
}

// This macro is a bit ugly, but trust me, it's better than the alternatives.
#define ROTATING_BODY_TAG_VALUE_CASE(value)                                    \
  case serialization::Frame::value:                                            \
    if (rotating_body_extension != nullptr) {                                  \
      return RotatingBody<                                                     \
                 Frame<Tag, serialization::Frame::value, true>>::              \
                 ReadFromMessage(*rotating_body_extension, parameters);        \
    } else {                                                                   \
      CHECK_NOTNULL(pre_brouwer_oblate_body_extension);                        \
      return OblateBody<                                                       \
                 Frame<Tag, serialization::Frame::value, true>>::              \
                 ReadFromMessage(                                              \
                     *pre_brouwer_oblate_body_extension, parameters);          \
    }

inline not_null<std::unique_ptr<MassiveBody>> MassiveBody::ReadFromMessage(
    serialization::MassiveBody const& message) {
  Parameters const parameters(GravitationalParameter::ReadFromMessage(
                                  message.gravitational_parameter()));

  // First see if we have an extension that has a frame and if so read the
  // frame.  Need to take care of pre-Brouwer compatibility.
  const google::protobuf::EnumValueDescriptor* enum_value_descriptor = nullptr;
  bool is_inertial = false;
  serialization::RotatingBody const* rotating_body_extension = nullptr;
  serialization::PreBrouwerOblateBody const* pre_brouwer_oblate_body_extension =
      nullptr;
  if (message.HasExtension(serialization::RotatingBody::rotating_body)) {
    rotating_body_extension =
        &message.GetExtension(serialization::RotatingBody::rotating_body);
    ReadFrameFromMessage(rotating_body_extension->frame(),
                         &enum_value_descriptor,
                         &is_inertial);
    CHECK(is_inertial);
  }
  if (message.HasExtension(
          serialization::PreBrouwerOblateBody::pre_brouwer_oblate_body)) {
    pre_brouwer_oblate_body_extension =
        &message.GetExtension(
            serialization::PreBrouwerOblateBody::pre_brouwer_oblate_body);
    ReadFrameFromMessage(pre_brouwer_oblate_body_extension->frame(),
                         &enum_value_descriptor,
                         &is_inertial);
    CHECK(is_inertial);
  }

  if (rotating_body_extension != nullptr ||
      pre_brouwer_oblate_body_extension != nullptr) {
    const google::protobuf::EnumDescriptor* enum_descriptor =
        enum_value_descriptor->type();
    {
      using Tag = serialization::Frame::PluginTag;
      if (enum_descriptor == google::protobuf::GetEnumDescriptor<Tag>()) {
        switch (static_cast<Tag>(enum_value_descriptor->number())) {
          ROTATING_BODY_TAG_VALUE_CASE(ALICE_SUN);
          ROTATING_BODY_TAG_VALUE_CASE(ALICE_WORLD);
          ROTATING_BODY_TAG_VALUE_CASE(BARYCENTRIC);
          ROTATING_BODY_TAG_VALUE_CASE(PRE_BOREL_BARYCENTRIC);
          ROTATING_BODY_TAG_VALUE_CASE(RENDERING);
          ROTATING_BODY_TAG_VALUE_CASE(WORLD);
          ROTATING_BODY_TAG_VALUE_CASE(WORLD_SUN);
        }
      }
    }
    {
      using Tag = serialization::Frame::SolarSystemTag;
      if (enum_descriptor == google::protobuf::GetEnumDescriptor<Tag>()) {
        switch (static_cast<Tag>(enum_value_descriptor->number())) {
          ROTATING_BODY_TAG_VALUE_CASE(ICRF_J2000_ECLIPTIC);
          ROTATING_BODY_TAG_VALUE_CASE(ICRF_J2000_EQUATOR);
        }
      }
    }
    {
      using Tag = serialization::Frame::TestTag;
      if (enum_descriptor == google::protobuf::GetEnumDescriptor<Tag>()) {
        switch (static_cast<Tag>(enum_value_descriptor->number())) {
          ROTATING_BODY_TAG_VALUE_CASE(TEST);
          ROTATING_BODY_TAG_VALUE_CASE(TEST1);
          ROTATING_BODY_TAG_VALUE_CASE(TEST2);
          ROTATING_BODY_TAG_VALUE_CASE(TEST3);
          ROTATING_BODY_TAG_VALUE_CASE(FROM);
          ROTATING_BODY_TAG_VALUE_CASE(THROUGH);
          ROTATING_BODY_TAG_VALUE_CASE(TO);
        }
      }
    }
    LOG(FATAL) << enum_descriptor->name();
    base::noreturn();
  } else {
    return std::make_unique<MassiveBody>(parameters);
  }
}

#undef ROTATING_BODY_TAG_VALUE_CASE

}  // namespace physics
}  // namespace principia
