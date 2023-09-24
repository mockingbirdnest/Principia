#pragma once

#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;

// Do not export these declarations, they are only exposed for tests.
template<typename Angle>
constexpr Angle one_π;

template<typename Angle>
constexpr Angle two_π;

template<double fractional_part_lower_bound,
         double fractional_part_upper_bound,
         typename Angle>
void ReduceAngle(Angle const& θ,
                 Angle& fractional_part,
                 std::int64_t& integer_part);

template<double fractional_part_lower_bound,
         double fractional_part_upper_bound,
         typename Angle>
Angle ReduceAngle(Angle const& θ);

}  // namespace internal

using internal::ReduceAngle;

}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia

#include "numerics/angle_reduction_body.hpp"
