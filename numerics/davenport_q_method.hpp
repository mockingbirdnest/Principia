#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/rotation.hpp"
#include "quantities/concepts.hpp"

namespace principia {
namespace numerics {
namespace _davenport_q_method {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_rotation;
using namespace principia::quantities::_concepts;

template<typename FromFrame, typename ToFrame, quantity Weight>
Rotation<FromFrame, ToFrame> DavenportQMethod(
    std::vector<Vector<double, FromFrame>> const& a,
    std::vector<Vector<double, ToFrame>> const& b,
    std::vector<Weight> const& weights);

}  // namespace internal

using internal::DavenportQMethod;

}  // namespace _davenport_q_method
}  // namespace numerics
}  // namespace principia

#include "numerics/davenport_q_method_body.hpp"
