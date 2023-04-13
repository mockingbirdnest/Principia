#pragma once

namespace principia {
namespace numerics {
namespace _ulp_distance {
namespace internal {

std::int64_t ULPDistance(double x, double y);

}  // namespace internal

using internal::ULPDistance;

}  // namespace _ulp_distance
}  // namespace numerics
}  // namespace principia

#include "numerics/ulp_distance_body.hpp"
