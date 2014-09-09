#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {

typedef Point<quantities::Time> Instant;
template<typename Frame>
using Displacement = Vector<quantities::Length, Frame>;
template<typename Frame>
using Position = Point<Displacement<Frame>>;
template<typename Frame>
using Velocity = Vector<quantities::Speed, Frame>;

}  // namespace geometry
}  // namespace principia
