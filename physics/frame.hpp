#pragma once

namespace principia {
namespace physics {

// The tag is just a way to unambiguously distinguish frames.
template<int tag>
class Frame {
  static int tag() { return tag; }
};

typedef Frame<0> InertialFrame;

}  // namespace physics
}  // namespace principia
