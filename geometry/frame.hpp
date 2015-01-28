#pragma once

#include "geometry/named_quantities.hpp"
#include "serialization/frame.pb.h"

namespace principia {
namespace geometry {

template<typename Tag, Tag tag, bool frame_is_inertial>
class Frame {
 public:
  static Position<Frame> const origin;
  static bool const is_inertial = frame_is_inertial;

  Frame() = delete;

  static void WriteToMessage(not_null<serialization::Frame*> const message);
  static void ReadFromMessage(serialization::Frame const& message);
};

// This frame should be used for objects whose reference frame cannot be known
// at compile time.
using UnknownInertialFrame = Frame<serialization::Frame::UnknownTag,
                                   serialization::Frame::UNKNOWN,
                                   true /*frame_is_inertial*/>;

}  // namespace geometry
}  // namespace principia

#include "geometry/frame_body.hpp"
