
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
         FrameTag tag_,
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

  static void WriteToMessage(not_null<serialization::Frame*> message);

  // Checks that the |message| matches the current type.
  static void ReadFromMessage(serialization::Frame const& message);
};

// The specializations below are for non-serializable frames.  They should be
// used as follows:
//   static int my_tag;
//   using MyFrame = Frame<void*, &my_tag>;

template<void* tag_,
         Inertia inertia_,
         Handedness handedness_>
struct Frame<void*, tag_, inertia_, handedness_> : not_constructible {
  static constexpr bool is_inertial = inertia_ == Inertial;
  static constexpr Inertia inertia = inertia_;
  static constexpr Handedness handedness = handedness_;

  static constexpr Position<Frame> const origin{};
  static constexpr Velocity<Frame> const unmoving{};
};

template<void* tag_>
struct Frame<void*, tag_, NonInertial, Handedness::Right> : not_constructible {
  static constexpr bool is_inertial = false;
  static constexpr Inertia inertia = NonInertial;
  static constexpr Handedness handedness = Handedness::Right;

  static constexpr Position<Frame> const origin{};
  static constexpr Velocity<Frame> const unmoving{};
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
