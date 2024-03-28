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

// This class is analogous to the pair which is its base class, except that it
// exports properly-named selectors.  It is implicitly convertible in both
// directions, so clients can generally ignore the difference.  Note however
// that creating a DegreesOfFreedom involves a copy so clients might want to use
// the base type (probably declared as |auto|) when they don't need to access
// the members.
template<typename Frame>
using DegreesOfFreedom = Pair<Position<Frame>, Velocity<Frame>>;

// This class is analogous to the vector class underlying DegreesOfFreedom,
// except that it exports properly-named selectors.  The same comments as above
// apply.
template<typename Frame>
using RelativeDegreesOfFreedom = Pair<Displacement<Frame>, Velocity<Frame>>;

}  // namespace internal

using internal::DegreesOfFreedom;
using internal::RelativeDegreesOfFreedom;

}  // namespace _degrees_of_freedom
}  // namespace physics
}  // namespace principia
