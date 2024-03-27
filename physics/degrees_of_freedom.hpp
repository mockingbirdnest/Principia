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

// Reopen the base namespace to make RelativeDegreesOfFreedom mappable.
namespace base {
namespace _mappable {
namespace internal {

using namespace principia::geometry::_pair;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;

template<typename Functor, typename Frame>
struct Mappable<Functor, RelativeDegreesOfFreedom<Frame>>
    : not_constructible {
  using type = Pair<decltype(std::declval<Functor>()(
                        std::declval<Displacement<Frame>>())),
                    decltype(std::declval<Functor>()(
                        std::declval<Velocity<Frame>>()))>;

  static type Do(Functor const& functor,
                 RelativeDegreesOfFreedom<Frame> const& relative);
};

}  // namespace internal
}  // namespace _mappable
}  // namespace base
}  // namespace principia

#include "physics/degrees_of_freedom_body.hpp"
