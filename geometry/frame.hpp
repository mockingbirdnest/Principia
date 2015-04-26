#pragma once

#include <cstdint>
#include <string>

#include "geometry/named_quantities.hpp"
#include "google/protobuf/descriptor.h"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

template<typename Tag_, Tag_ tag, bool is_inertial>
class Frame {
 public:
  using Tag = Tag_;
  static Position<Frame> const origin;
  static Tag const tag = tag;
  static bool const is_inertial = is_inertial;

  Frame() = delete;

  static void WriteToMessage(not_null<serialization::Frame*> const message);

  // Checks that the |message| matches the current type.
  static void ReadFromMessage(serialization::Frame const& message);
};

// Extracts enough information from the |message| to contruct a |Frame| type.
void ReadFrameFromMessage(
    serialization::Frame const& message,
    not_null<google::protobuf::EnumValueDescriptor const**> const
        enum_value_descriptor,
    not_null<bool*> const is_inertial);

}  // namespace geometry
}  // namespace principia

#include "geometry/frame_body.hpp"
