#pragma once

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Vector;
using principia::quantities::Length;
using principia::quantities::Speed;

namespace principia {
namespace physics {

template<typename Frame>
struct DegreesOfFreedom {
  DegreesOfFreedom(Vector<Length, Frame> const& position,
                   Vector<Speed, Frame> const& velocity);
  Vector<Length, Frame> const position;
  Vector<Speed, Frame> const velocity;
};

}  // namespace physics
}  // namespace principia

#include "physics/degrees_of_freedom_body.hpp"
