#pragma once

#include "numerics/finite_difference.hpp"

#include "numerics/finite_difference.mathematica.h"

namespace principia {
namespace numerics {
namespace internal_finite_difference {

template<typename Value, typename Argument, std::size_t n>
Derivative<Value, Argument> FiniteDifference(
    std::array<Value, n> const& values,
    Argument const& step,
    int const offset) {
  double const* const numerators = std::get<n - 1>(Numerators)[offset];
  constexpr double denominator = Denominators[n - 1];
  Difference<Value> sum{};
  if (n % 2 == 1 && offset == (n - 1) / 2) {
    // For the central difference formula, aᵢ = - aₙ₋ᵢ₋₁; in particular, for
    // i = (n - 1) / 2 (the central coefficient), aᵢ = -aᵢ: the central value is
    // unused.
    // We thus evaluate the sum Σᵢ aᵢ f(xᵢ), with i runnning from 0 to n - 1, as
    // Σⱼ aⱼ (f(xⱼ) - f(xₙ₋ⱼ₋₁)), with j running from 0 to (n - 3) / 2.  Which
    // we cannot write naively because n is unsigned.
    for (int j = 0; 2 * j + 3 <= n; ++j) {
      sum += numerators[j] * (values[j] - values[n - j - 1]);
    }
    return sum / (denominator * step);
  } else {
    // In the general case, we evaluate the sum Σᵢ aᵢ f(xᵢ), with Σᵢ aᵢ = 0,
    // where the sums over i run from 0 to n - 1, as
    //   Σⱼ (Σₖ aₖ) (f(xⱼ) - f(xⱼ₊₁)),
    // where the sum over j runs from 0 to n - 2, and the sum over
    // k runs from 0 to j.
    double numerator = 0;
    for (int j = 0; j + 2 <= n; ++j) {
      numerator += numerators[j];
      sum += numerator * (values[j] - values[j + 1]);
    }
  }
  return sum / (denominator * step);
}

}  // namespace internal_finite_difference

using internal_finite_difference::FiniteDifference;

}  // namespace numerics
}  // namespace principia
