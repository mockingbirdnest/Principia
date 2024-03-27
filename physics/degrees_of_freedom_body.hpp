#pragma once

#include "physics/degrees_of_freedom.hpp"

#include <string>
#include "geometry/barycentre_calculator.hpp"
#include "geometry/pair.hpp"

namespace principia {
namespace physics {
namespace _pair {
namespace internal {

using namespace _degrees_of_freedom;

// REMOVE BEFORE FLIGHT Move these up into Pair.
template<typename Frame>
std::string DebugString(DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  return "{" + DebugString(degrees_of_freedom.position()) + ", " +
         DebugString(degrees_of_freedom.velocity()) + "}";
}

template<typename Frame>
std::string DebugString(
    RelativeDegreesOfFreedom<Frame> const& relative_degrees_of_freedom) {
  return "{" + DebugString(relative_degrees_of_freedom.displacement()) + ", " +
         DebugString(relative_degrees_of_freedom.velocity()) + "}";
}

template<typename Frame>
std::ostream& operator<<(std::ostream& out,
                         DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  out << DebugString(degrees_of_freedom);
  return out;
}

template<typename Frame>
std::ostream& operator<<(
    std::ostream& out,
    RelativeDegreesOfFreedom<Frame> const& relative_degrees_of_freedom) {
  out << DebugString(relative_degrees_of_freedom);
  return out;
}

}  // namespace internal
}  // namespace _pair
}  // namespace physics

namespace base {
namespace _mappable {
namespace internal {

using namespace principia::geometry::_pair;
using namespace principia::physics::_degrees_of_freedom;

template<typename Functor, typename Frame>
typename Mappable<Functor, RelativeDegreesOfFreedom<Frame>>::type
Mappable<Functor, RelativeDegreesOfFreedom<Frame>>::Do(
    Functor const& functor,
    RelativeDegreesOfFreedom<Frame> const& relative) {
  return type(functor(relative.displacement()), functor(relative.velocity()));
}

}  // namespace internal
}  // namespace _mappable
}  // namespace base
}  // namespace principia
