#pragma once

#include <cstdint>
#include <string>

#include "base/not_constructible.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/descriptor.h"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _frame {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;

// The enumerators of |FrameMotion| are ordered from most restrictive to least
// restrictive; m1 <= m2 means that m1 satisfies the requirements of m2.
enum FrameMotion {
  Inertial,
  NonRotating,
  Arbitrary,
};

enum class Handedness {
  Left,
  Right,
};

// The frame is serializable if and only if FrameTag is a protocol buffer enum.
// To declare a local frame that does not need serialization, use the following
// pattern:
//   using MyFrame = Frame<struct MyFrameTag>;
//
// or:
//   using MyFrame = Frame<struct MyFrameTag, Inertial>;
//
// By default, the frame is arbitrary and right-handed.
//
// A non-serializable frame misses the ReadFromMessage method but has a
// WriteToMessage method which fails at execution.  The reason is that
// WriteToMessage needs to be virtual for some classes, and that doesn't
// interact well with SFINAE.  Since it's uncommon to write without reading,
// that should be sufficient to detect at compile-time attempts at serializing a
// non-serializable frame.
template<typename FrameTag,
         FrameMotion motion_ = Arbitrary,
         Handedness handedness_ = Handedness::Right,
         std::conditional_t<std::is_enum_v<FrameTag>,
                            FrameTag, std::nullptr_t> tag_ =
             std::conditional_t<std::is_enum_v<FrameTag>,
                                FrameTag, std::nullptr_t>{}>
struct Frame : not_constructible {
  static constexpr bool is_inertial = motion_ == Inertial;
  static constexpr bool may_rotate = motion_ == Arbitrary;
  static constexpr FrameMotion motion = motion_;
  static constexpr Handedness handedness = handedness_;

  static const Position<Frame> origin;
  static const Velocity<Frame> unmoving;
  static const AngularVelocity<Frame> nonrotating;

  using Tag = FrameTag;
  static constexpr bool has_tag = std::is_enum_v<FrameTag>;
  template<typename T = FrameTag,
           typename = std::enable_if_t<std::is_enum_v<T>>>
  static constexpr FrameTag tag();

  static void WriteToMessage(not_null<serialization::Frame*> message);

  // Checks that the |message| matches the current type.
  template<
      typename T = FrameTag,
      typename = std::enable_if_t<google::protobuf::is_proto_enum<T>::value>>
  static void ReadFromMessage(serialization::Frame const& message);
};

// Extracts enough information from the |message| to contruct a |Frame| type.
void ReadFrameFromMessage(
    serialization::Frame const& message,
    google::protobuf::EnumValueDescriptor const*& enum_value_descriptor,
    bool& is_inertial);

}  // namespace internal

using internal::Arbitrary;
using internal::Frame;
using internal::Handedness;
using internal::Inertial;
using internal::NonRotating;
using internal::ReadFrameFromMessage;

}  // namespace _frame
}  // namespace geometry
}  // namespace principia

namespace principia::geometry {
using namespace principia::geometry::_frame;
}  // namespace principia::geometry

#include "geometry/frame_body.hpp"
