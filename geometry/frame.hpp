
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

enum FrameMotion {
  Inertial,
  NonInertial,
};

enum class Handedness {
  Left,
  Right,
};

// The frame is serializable if and only if FrameTag is a protocol buffer enum.
// To declare a local frame that does not need serialization, use the following
// pattern:
//   using MyFrame = Frame<enum class MyFrameTag>;
//
// or:
//   using MyFrame = Frame<enum class MyFrameTag, Inertial>;
//
// By default, the frame is non-inertial and right-handed.
//
// A non-serializable frame misses the ReadFromMessage method but has a
// WriteToMessage method which fails at execution.  The reason is that
// WriteToMessage needs to be virtual for some classes, and that doesn't
// interact well with SFINAE.  Since it's uncommon to write without reading,
// that should be sufficient to detect at compile-time attempts at serializing a
// non-serializable frame.
template<typename FrameTag,
         FrameMotion motion_ = NonInertial,
         Handedness handedness_ = Handedness::Right,
         FrameTag tag_ = FrameTag{}>
struct Frame : not_constructible {
  static constexpr bool is_inertial = motion_ == Inertial;
  static constexpr FrameMotion motion = motion_;
  static constexpr Handedness handedness = handedness_;

  static const Position<Frame> origin;
  static const Velocity<Frame> unmoving;

  using Tag = FrameTag;
  static constexpr Tag tag = tag_;

  static constexpr bool is_serializable =
      google::protobuf::is_proto_enum<FrameTag>::value;

  static void WriteToMessage(not_null<serialization::Frame*> message);

  // Checks that the |message| matches the current type.
  template<typename = std::enable_if_t<is_serializable>>
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
