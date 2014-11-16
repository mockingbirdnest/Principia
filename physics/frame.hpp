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

// A helper for declaring tags without conflicts.
// TODO(phl): Ideally we'd like a constexpr function hashing the
// (__FILE__, __LINE__) pair.
#define UNIQUE_TAG (1000 * __COUNTER__ + __LINE__)

}  // namespace physics
}  // namespace principia

#include "physics/frame_body.hpp"
