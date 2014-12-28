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

}  // namespace physics
}  // namespace principia

#include "physics/frame_body.hpp"
