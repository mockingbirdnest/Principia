#pragma once

#include "numerics/finite_difference.hpp"

#include "numerics/double_precision.hpp"
#include "numerics/finite_difference.mathematica.h"

namespace principia {
namespace numerics {
namespace internal_finite_difference {

template<typename Value, typename Argument, int n>
Derivative<Value, Argument> FiniteDifference(
    FixedVector<Value, n> const& values,
    Argument const& step,
    int offset) {
  double const* const numerators = std::get<n - 1>(Numerators)[offset];
  constexpr double denominator = Denominators[n - 1];
  if (n % 2 == 1 && offset == (n - 1) / 2) {
    // For the central difference formula, aᵢ = - aₙ₋ᵢ₋₁; in particular, for
    // i = (n - 1) / 2 (the central coefficient), aᵢ = -aᵢ: the central value is
    // unused.
    // We thus evaluate the sum Σᵢ aᵢ f(xᵢ), with i runnning from 0 to n - 1, as
    // Σⱼ aⱼ (f(xⱼ) - f(xₙ₋ⱼ₋₁)), with j running from 0 to (n - 3) / 2.
    DoublePrecision<Difference<Value>> sum{};
    for (int j = 0; j <= (n - 3) / 2; ++j) {
      sum += TwoProduct(numerators[j], values[j] - values[n - j - 1]);
    }
    return sum.value / (denominator * step);
  } else {
    // In the general case, we evaluate the sum Σᵢ aᵢ f(xᵢ), with Σᵢ aᵢ = 0,
    // where the sums over i run from 0 to n - 1, as
    //   Σⱼ (Σₖ aₖ) (f(xⱼ) - f(xⱼ₊₁)),
    // where the sum over j runs from 0 to n - 2, and the sum over
    // k runs from 0 to j.
    double numerator = 0;
    DoublePrecision<Difference<Value>> sum{};
    for (int j = 0; j <= n - 2; ++j) {
      numerator += numerators[j];
      Difference<Value> difference = values[j] - values[j + 1];
      sum += TwoProduct(numerator, difference);
    }
    return sum.value / (denominator * step);
  }
}

}  // namespace internal_finite_difference

using internal_finite_difference::FiniteDifference;

}  // namespace numerics
}  // namespace principia
