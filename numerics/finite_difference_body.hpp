#pragma once

#include "numerics/finite_difference.hpp"

#include "numerics/double_precision.hpp"

namespace principia {
namespace numerics {
namespace internal_finite_difference {

// Represents a linear combination of values with coefficients |numerators[i]|
// divided by the common |denominator|.
// The sum of the entries in |numerators| must be 0.

constexpr int max_order = 5;

constexpr FixedLowerTriangularMatrix<double, max_order + 1> forward_differences{
    {{
        0,
        -1, 1,
        -3.0 / 2,        2,        -1.0 / 2,
        -11.0/6,   3.0, -3.0/2.0, 1.0/3.0,
         -25.0/12,  4,   -3, 4.0/3, -1.0/4,
         -137, 300, -300, 200, -75, 12,
    }}};

constexpr double denominator = 60;

static int IMPL;

template<typename Value, typename Argument, int n>
Derivative<Value, Argument> FiniteDifference(
    FixedVector<Value, n> const& values,
    Argument const& step,
    int offset) {
  CHECK_EQ(offset, 0);
  constexpr double const* numerators = forward_differences[n - 1];
  LOG(ERROR)<<IMPL;

if (IMPL == 1) {
  Difference<Value> sum{};
  for (int j = 0; j <= n-1; ++j) {
    sum += (numerators[j] / denominator) * values[j];
  }
  return sum / (1 * step);
} else if  (IMPL == 2) {
  // We evaluate the sum Σᵢ aᵢ f(xᵢ), with Σᵢ aᵢ = 0, where the sums over i runs
  // from 0 to n - 1, as
  //   Σⱼ (Σₖ aₖ) (f(xⱼ) - f (xⱼ₊₁)),
  // where the sum over j runs from 0 to n - 2, and the sum over
  // k runs from 0 to j.
  double numerator = 0;
  Difference<Value> sum{};
  for (int j = 0; j <= n - 2; ++j) {
    numerator += numerators[j] / denominator;
    Difference<Value> difference = values[j] - values[j + 1];
    sum += numerator * difference;
  }
  return sum / (1 * step);
} else if  (IMPL == 3) {
  // We evaluate the sum Σᵢ aᵢ f(xᵢ), with Σᵢ aᵢ = 0, where the sums over i runs
  // from 0 to n - 1, as
  //   Σⱼ (Σₖ aₖ) (f(xⱼ) - f (xⱼ₊₁)),
  // where the sum over j runs from 0 to n - 2, and the sum over
  // k runs from 0 to j.
  double numerator = 0;
  DoublePrecision<Difference<Value>> sum{};
  for (int j = 0; j <= n - 2; ++j) {
    numerator += numerators[j] / denominator;
    Difference<Value> difference = values[j] - values[j + 1];
    sum.Increment(numerator * difference);
  }
  return sum.value / step;
} else { CHECK_EQ(IMPL, 4);
  // We evaluate the sum Σᵢ aᵢ f(xᵢ), with Σᵢ aᵢ = 0, where the sums over i runs
  // from 0 to n - 1, as
  //   Σⱼ (Σₖ aₖ) (f(xⱼ) - f (xⱼ₊₁)),
  // where the sum over j runs from 0 to n - 2, and the sum over
  // k runs from 0 to j.
  double numerator = 0;
  DoublePrecision<Difference<Value>> sum{};
  for (int j = 0; j <= n - 2; ++j) {
    numerator += numerators[j];
    Difference<Value> difference = values[j] - values[j + 1];
    sum.Increment(numerator * difference);
  }
  return sum.value / (denominator * step);
}
}

}  // namespace internal_finite_difference

using internal_finite_difference::FiniteDifference;

}  // namespace numerics
}  // namespace principia
