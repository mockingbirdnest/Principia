#pragma once

#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;

// Clang up to version 17 does not support non-type template parameters that are
// floats.  So we wrap the float in a silly struct instead.
#if (PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL) && \
    __clang_major__ <= 17
struct DoubleWrapper {
  constexpr DoubleWrapper(double const d) : d(d) {}  // NOLINT(runtime/explicit)
  constexpr DoubleWrapper(int const i) : d(i) {}  // NOLINT(runtime/explicit)
  double d;
};
#else
using DoubleWrapper = double;
#endif

// Do not export these declarations, they are only exposed for tests.
template<typename Angle>
constexpr Angle one_π;

template<typename Angle>
constexpr Angle two_π;

// Performs Payne-Hanek reduction of `x`.  `x_reduced` is in [-π / 4, π / 4],
// and the `quadrant` ranges from 0 to 3.  0 is a quadrant centered on the
// positive x axis, 1 on the positive y axis, and so on.
template<std::int64_t precision>
void PayneHanek(Angle const& x,
                DoublePrecision<Angle>& x_reduced,
                std::int64_t& quadrant);

// If [fractional_part_lower_bound, fractional_part_upper_bound] covers 2π or
// more, the reduction is modulo 2π.  If it covers only π, the reduction is
// modulo π.

// Return false if the argument is too large.
template<DoubleWrapper fractional_part_lower_bound,
         DoubleWrapper fractional_part_upper_bound,
         typename Angle>
bool ReduceAngle(Angle const& θ,
                 Angle& fractional_part,
                 std::int64_t& integer_part);

// Fails if the argument is too large.
template<DoubleWrapper fractional_part_lower_bound,
         DoubleWrapper fractional_part_upper_bound,
         typename Angle>
Angle ReduceAngle(Angle const& θ);

}  // namespace internal

using internal::PayneHanek;
using internal::ReduceAngle;

}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia

#include "numerics/angle_reduction_body.hpp"
