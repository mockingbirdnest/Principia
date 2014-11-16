#pragma once

#include "geometry/named_quantities.hpp"

using principia::geometry::Position;

namespace principia {
namespace physics {

template<int tag, bool is_inertial>
class Frame {
 public:
  static Position<Frame> const origin;
  static bool const is_inertial = is_inertial;

  Frame() = delete;
};

}  // namespace physics
}  // namespace principia

#include "physics/frame_body.hpp"