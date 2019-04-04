#pragma once

#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_finite_difference {

using quantities::Derivative;
using quantities::Difference;

// Given n equally-spaced values f(xᵢ) = f(x₀ + i h),
// approximates the derivative f′(xⱼ) at xⱼ for j = |offset|, where 0 ≤ j < n.
// Special values of |offset| are:
// — |offset = 0|: forward difference;
// — |offset = (n - 1) / 2|, for odd n:
//   central difference—in this case, the middle value is unused;
// — |offset = n - 1|: backward difference.
// The error on the derivative is 𝒪(hⁿ⁻¹) as h → 0.
// If f is a polynomial of degree less than or equal to n - 1, the result is
// exact up to rounding errors.
template<typename Value, typename Argument, int n>
Derivative<Value, Argument> FiniteDifference(
    FixedVector<Value, n> const& values,
    Argument const& step,
    int offset);

}  // namespace internal_finite_difference

using internal_finite_difference::FiniteDifference;

}  // namespace numerics
}  // namespace principia

#include "numerics/finite_difference_body.hpp"
