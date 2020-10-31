
#pragma once

#include "physics/massive_body.hpp"

#include <string>
#include <utility>

#include "base/macros.hpp"
#include "geometry/frame.hpp"
#include "glog/logging.h"
#include "google/protobuf/descriptor.h"
#include "physics/oblate_body.hpp"
#include "quantities/constants.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace physics {
namespace internal_massive_body {

using quantities::constants::GravitationalConstant;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::ReadFrameFromMessage;

inline MassiveBody::Parameters::Parameters(
    GravitationalParameter const& gravitational_parameter)
    : Parameters(/*name=*/"", gravitational_parameter) {}

inline MassiveBody::Parameters::Parameters(
    std::string name,
    GravitationalParameter const& gravitational_parameter)
    : name_(std::move(name)),
      gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant) {
  CHECK_NE(gravitational_parameter, GravitationalParameter())
      << "Massive body cannot have zero gravitational parameter";
}

inline MassiveBody::Parameters::Parameters(Mass const& mass)
    : Parameters(/*name=*/"", mass) {}

inline MassiveBody::Parameters::Parameters(std::string name,
                                           Mass const& mass)
    : name_(std::move(name)),
      gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass) {
  CHECK_NE(mass, Mass()) << "Massive body cannot have zero mass";
}

inline GravitationalParameter const&
MassiveBody::Parameters::gravitational_parameter() const {
  return gravitational_parameter_;
}

inline MassiveBody::MassiveBody(Parameters parameters)
    : parameters_(std::move(parameters)) {}

inline std::string const& MassiveBody::name() const {
  return parameters_.name_;
}

inline GravitationalParameter const&
MassiveBody::gravitational_parameter() const {
  return parameters_.gravitational_parameter_;
}

inline Mass const& MassiveBody::mass() const {
  return parameters_.mass_;
}

inline Length MassiveBody::min_radius() const {
  return Length();
}

inline Length MassiveBody::mean_radius() const {
  return Length();
}

inline Length MassiveBody::max_radius() const {
  return Length();
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
  message->set_name(parameters_.name_);
  parameters_.gravitational_parameter_.WriteToMessage(
      message->mutable_gravitational_parameter());
}

inline not_null<std::unique_ptr<MassiveBody>> MassiveBody::ReadFromMessage(
    serialization::Body const& message) {
  CHECK(message.has_massive_body());
  return ReadFromMessage(message.massive_body());
}

// This macro is a bit ugly, but trust me, it's better than the alternatives.
#define ROTATING_BODY_TAG_VALUE_CASE(value)                   \
  case serialization::Frame::value:                           \
    CHECK_NOTNULL(rotating_body_extension);                   \
    return RotatingBody<Frame<Tag,                            \
                              Inertial,                       \
                              Handedness::Right,              \
                              serialization::Frame::value>>:: \
        ReadFromMessage(*rotating_body_extension, parameters)

inline not_null<std::unique_ptr<MassiveBody>> MassiveBody::ReadFromMessage(
    serialization::MassiveBody const& message) {
  Parameters const parameters(message.name(),
                              GravitationalParameter::ReadFromMessage(
                                  message.gravitational_parameter()));

  // First see if we have an extension that has a frame and if so read the
  // frame.
  const google::protobuf::EnumValueDescriptor* enum_value_descriptor = nullptr;
  bool is_inertial = false;
  serialization::RotatingBody const* rotating_body_extension = nullptr;
  if (message.HasExtension(serialization::RotatingBody::extension)) {
    rotating_body_extension =
        &message.GetExtension(serialization::RotatingBody::extension);
    ReadFrameFromMessage(rotating_body_extension->frame(),
                         enum_value_descriptor,
                         is_inertial);
    CHECK(is_inertial);
  }

  if (rotating_body_extension != nullptr) {
    const google::protobuf::EnumDescriptor* enum_descriptor =
        enum_value_descriptor->type();
    {
      using Tag = serialization::Frame::PhysicsTag;
      if (enum_descriptor == google::protobuf::GetEnumDescriptor<Tag>()) {
        switch (static_cast<Tag>(enum_value_descriptor->number())) {
          ROTATING_BODY_TAG_VALUE_CASE(FRENET);
          ROTATING_BODY_TAG_VALUE_CASE(PRINCIPAL_AXES);
        }
      }
    }
    {
      using Tag = serialization::Frame::PluginTag;
      if (enum_descriptor == google::protobuf::GetEnumDescriptor<Tag>()) {
        switch (static_cast<Tag>(enum_value_descriptor->number())) {
          ROTATING_BODY_TAG_VALUE_CASE(ALICE_SUN);
          ROTATING_BODY_TAG_VALUE_CASE(ALICE_WORLD);
          ROTATING_BODY_TAG_VALUE_CASE(APPARENT);
          ROTATING_BODY_TAG_VALUE_CASE(APPARENT_WORLD);
          ROTATING_BODY_TAG_VALUE_CASE(BARYCENTRIC);
          ROTATING_BODY_TAG_VALUE_CASE(BODY_WORLD);
          ROTATING_BODY_TAG_VALUE_CASE(CAMERA);
          ROTATING_BODY_TAG_VALUE_CASE(CAMERA_REFERENCE);
          ROTATING_BODY_TAG_VALUE_CASE(CELESTIAL_SPHERE);
          ROTATING_BODY_TAG_VALUE_CASE(ECCENTRIC_PART);
          ROTATING_BODY_TAG_VALUE_CASE(MAIN_BODY_CENTRED);
          ROTATING_BODY_TAG_VALUE_CASE(NAVBALL);
          ROTATING_BODY_TAG_VALUE_CASE(NAVIGATION);
          ROTATING_BODY_TAG_VALUE_CASE(NON_ROTATING_PILE_UP);
          ROTATING_BODY_TAG_VALUE_CASE(PILE_UP_PRINCIPAL_AXES);
          ROTATING_BODY_TAG_VALUE_CASE(RIGID_PART);
          ROTATING_BODY_TAG_VALUE_CASE(WORLD);
          ROTATING_BODY_TAG_VALUE_CASE(WORLD_SUN);
        }
      }
    }
    {
      using Tag = serialization::Frame::SolarSystemTag;
      if (enum_descriptor == google::protobuf::GetEnumDescriptor<Tag>()) {
        switch (static_cast<Tag>(enum_value_descriptor->number())) {
          ROTATING_BODY_TAG_VALUE_CASE(GCRS);
          ROTATING_BODY_TAG_VALUE_CASE(ICRS);
          ROTATING_BODY_TAG_VALUE_CASE(ITRS);
          ROTATING_BODY_TAG_VALUE_CASE(SKY);
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

}  // namespace internal_massive_body
}  // namespace physics
}  // namespace principia
