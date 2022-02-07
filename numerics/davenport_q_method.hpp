#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/quaternion.hpp"

namespace principia {
namespace numerics {
namespace internal_davenport_q_method {

using geometry::Quaternion;
using geometry::Vector;

template<typename Weight, typename Frame>
Quaternion DavenportQMethod(std::vector<Vector<double, Frame>> const& a,
                            std::vector<Vector<double, Frame>> const& b,
                            std::vector<Weight> const& weights);

}  // namespace internal_davenport_q_method

using internal_davenport_q_method::DavenportQMethod;

}  // namespace numerics
}  // namespace principia

#include "numerics/davenport_q_method_body.hpp"
