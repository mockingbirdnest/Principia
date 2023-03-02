#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/quaternion.hpp"

namespace principia {
namespace numerics {
namespace internal_davenport_q_method {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_rotation;

template<typename FromFrame, typename ToFrame, typename Weight>
Rotation<FromFrame, ToFrame> DavenportQMethod(
    std::vector<Vector<double, FromFrame>> const& a,
    std::vector<Vector<double, ToFrame>> const& b,
    std::vector<Weight> const& weights);

}  // namespace internal_davenport_q_method

using internal_davenport_q_method::DavenportQMethod;

}  // namespace numerics
}  // namespace principia

#include "numerics/davenport_q_method_body.hpp"
