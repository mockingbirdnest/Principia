#pragma once

#include "geometry/direct_sum.hpp"
#include "geometry/space.hpp"

namespace principia {
namespace physics {
namespace _degrees_of_freedom {
namespace internal {

using namespace principia::geometry::_direct_sum;
using namespace principia::geometry::_space;

template<typename Frame>
using DegreesOfFreedom = DirectSum<Position<Frame>, Velocity<Frame>>;

template<typename Frame>
using RelativeDegreesOfFreedom =
    DirectSum<Displacement<Frame>, Velocity<Frame>>;

}  // namespace internal

using internal::DegreesOfFreedom;
using internal::RelativeDegreesOfFreedom;

}  // namespace _degrees_of_freedom
}  // namespace physics
}  // namespace principia
