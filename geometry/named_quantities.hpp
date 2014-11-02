#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {

using Instant = Point<quantities::Time>;
template<typename Frame>
using Displacement = Vector<quantities::Length, Frame>;
template<typename Frame>
using Position = Point<Displacement<Frame>>;
template<typename Frame>
using Velocity = Vector<quantities::Speed, Frame>;

template<typename Frame>
using AngularVelocity = Bivector<quantities::AngularFrequency, Frame>;

}  // namespace geometry
}  // namespace principia
