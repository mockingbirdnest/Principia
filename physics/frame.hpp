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

template<int tag, bool is_inertial>
Position<Frame<tag, is_inertial>> const Frame<tag, is_inertial>::origin;

template<int tag, bool is_inertial>
bool const Frame<tag, is_inertial>::is_inertial;

}  // namespace physics
}  // namespace principia