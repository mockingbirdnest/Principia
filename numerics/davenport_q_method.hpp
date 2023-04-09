#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/rotation.hpp"

namespace principia {
namespace numerics {
namespace _davenport_q_method {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_rotation;

template<typename FromFrame, typename ToFrame, typename Weight>
Rotation<FromFrame, ToFrame> DavenportQMethod(
    std::vector<Vector<double, FromFrame>> const& a,
    std::vector<Vector<double, ToFrame>> const& b,
    std::vector<Weight> const& weights);

}  // namespace internal

using internal::DavenportQMethod;

}  // namespace _davenport_q_method
}  // namespace numerics
}  // namespace principia

namespace principia::numerics {
using namespace principia::numerics::_davenport_q_method;
}  // namespace principia::numerics

#include "numerics/davenport_q_method_body.hpp"
