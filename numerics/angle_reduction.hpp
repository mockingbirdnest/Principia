#pragma once

#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;

// Argument reduction: angle = fractional_part + integer_part * π where
// fractional_part is in [-π/2, π/2].
void Reduce(Angle const& angle,
            Angle& fractional_part,
            std::int64_t& integer_part);

DoublePrecision<Angle> Mod2π(DoublePrecision<Angle> const& θ);

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
}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia

#include "numerics/angle_reduction_body.hpp"
