#pragma once

#include "geometry/named_quantities.hpp"

using principia::geometry::Position;

namespace principia {
namespace physics {

template<typename Tag, Tag tag, bool frame_is_inertial>
class Frame {
 public:
  static Position<Frame> const origin;
  static bool const is_inertial = frame_is_inertial;

  Frame() = delete;
};

enum UncheckedTag {
  kUnchecked,
};

// This frame should be used for objects whose reference frame cannot be known
// at compile time.
using UncheckedInertialFrame =
    Frame<UncheckedTag, UncheckedTag::kUnchecked, true>;

}  // namespace physics
}  // namespace principia

#include "physics/frame_body.hpp"
