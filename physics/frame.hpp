#pragma once

#include "geometry/named_quantities.hpp"

using principia::geometry::Position;

namespace principia {
namespace physics {

template<typename Tag, Tag tag, bool is_inertialLOLZ>
class Frame {
 public:
  static Position<Frame> const origin;
  static bool const is_inertial = is_inertialLOLZ;

  Frame() = delete;
};

}  // namespace physics
}  // namespace principia

#include "physics/frame_body.hpp"
