#pragma once

#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;

// Performs Payne-Hanek reduction of `x`.  `x_reduced` is in [-π / 4, π / 4],
// and the `quadrant` ranges from 0 to 3.  0 is a quadrant centered on the
// positive x axis, 1 on the positive y axis, and so on.  This function is
// usable with types `Angle` and `double`.
template<std::int64_t precision, typename Angle>
void PayneHanek(Angle const& x,
                DoublePrecision<Angle>& x_reduced,
                std::int64_t& quadrant);

// If [fractional_part_lower_bound, fractional_part_upper_bound] covers 2π or
// more, the reduction is modulo 2π.  If it covers only π, the reduction is
// modulo π.  The `integer_part` is a double to account for very large angles.

template<double fractional_part_lower_bound,
         double fractional_part_upper_bound>
Angle ReduceAngle(Angle const& θ);

template<double fractional_part_lower_bound,
         double fractional_part_upper_bound>
void ReduceAngle(Angle const& θ,
                 Angle& fractional_part,
                 double& integer_part);

}  // namespace internal

using internal::PayneHanek;
using internal::ReduceAngle;

}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia

#include "numerics/angle_reduction_body.hpp"
