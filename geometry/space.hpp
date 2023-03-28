#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace _space {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Frame>
using Displacement = Vector<Length, Frame>;

template<typename Frame>
using Position = Point<Displacement<Frame>>;

template<typename Frame>
using Velocity = Vector<Speed, Frame>;

template<typename Frame>
using AngularVelocity = Bivector<AngularFrequency, Frame>;

}  // namespace internal

using internal::AngularVelocity;
using internal::Displacement;
using internal::Position;
using internal::Velocity;

}  // namespace _space
}  // namespace geometry
}  // namespace principia
