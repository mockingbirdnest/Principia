#pragma once

#include <array>

#include "quantities/arithmetic.hpp"

namespace principia {
namespace numerics {
namespace _finite_difference {
namespace internal {

using namespace principia::quantities::_arithmetic;

// Given n equally-spaced values f(xᵢ) = f(x₀ + i h),
// approximates the derivative f′(xⱼ) at xⱼ for j = `offset`, where 0 ≤ j < n.
// Special values of `offset` are:
// — `offset = 0`: forward difference;
// — `offset = (n - 1) / 2`, for odd n:
//   central difference—in this case, the middle value is unused;
// — `offset = n - 1`: backward difference.
// The error on the derivative is 𝒪(hⁿ⁻¹) as h → 0.
// If f is a polynomial of degree less than or equal to n - 1, the result is
// exact up to rounding errors.
template<typename Value, typename Argument, std::size_t n>
Derivative<Value, Argument> FiniteDifference(
    std::array<Value, n> const& values,
    Argument const& step,
    int offset);

}  // namespace internal

using internal::FiniteDifference;

}  // namespace _finite_difference
}  // namespace numerics
}  // namespace principia

#include "numerics/finite_difference_body.hpp"
