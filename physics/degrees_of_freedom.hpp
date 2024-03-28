#pragma once

#include <string>
#include <vector>

#include "base/concepts.hpp"
#include "base/not_constructible.hpp"
#include "geometry/pair.hpp"
#include "geometry/space.hpp"

namespace principia {
namespace physics {
namespace _degrees_of_freedom {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_not_constructible;
using namespace principia::geometry::_pair;
using namespace principia::geometry::_space;

template<typename Frame>
using DegreesOfFreedom = Pair<Position<Frame>, Velocity<Frame>>;

template<typename Frame>
using RelativeDegreesOfFreedom = Pair<Displacement<Frame>, Velocity<Frame>>;

}  // namespace internal

using internal::DegreesOfFreedom;
using internal::RelativeDegreesOfFreedom;

}  // namespace _degrees_of_freedom
}  // namespace physics
}  // namespace principia
