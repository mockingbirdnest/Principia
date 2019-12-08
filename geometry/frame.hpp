
#pragma once

#include <cstdint>
#include <string>

#include "base/not_constructible.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/descriptor.h"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_frame {

using base::not_constructible;
using base::not_null;

enum Inertia {
  Inertial,
  NonInertial,
};

enum class Handedness {
  Left,
  Right,
};

template<typename FrameTag,
         FrameTag tag_ = FrameTag{},
         Inertia inertia_ = NonInertial,
         Handedness handedness_ = Handedness::Right>
struct Frame : not_constructible {
  static constexpr bool is_inertial = inertia_ == Inertial;
  static constexpr Inertia inertia = inertia_;
  static constexpr Handedness handedness = handedness_;

  static constexpr Position<Frame> const origin{};
  static constexpr Velocity<Frame> const unmoving{};

  using Tag = FrameTag;
  static constexpr Tag tag = tag_;

  template<typename = std::enable_if<
               google::protobuf::is_proto_enum<FrameTag>::value>>
  static void WriteToMessage(not_null<serialization::Frame*> message);

  // Checks that the |message| matches the current type.
  template<typename = std::enable_if<
               google::protobuf::is_proto_enum<FrameTag>::value>>
  static void ReadFromMessage(serialization::Frame const& message);
};

// Extracts enough information from the |message| to contruct a |Frame| type.
void ReadFrameFromMessage(
    serialization::Frame const& message,
    google::protobuf::EnumValueDescriptor const*& enum_value_descriptor,
    bool& is_inertial);

}  // namespace internal_frame

using internal_frame::Frame;
using internal_frame::Handedness;
using internal_frame::Inertial;
using internal_frame::NonInertial;
using internal_frame::ReadFrameFromMessage;

}  // namespace geometry
}  // namespace principia

#include "geometry/frame_body.hpp"
