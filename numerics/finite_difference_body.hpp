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
  CHECK_EQ(offset, 0);
  constexpr FixedMatrix<double, n, n> numerators = std::get<n - 1>(Numerators);
  constexpr double denominator = denominators[n - 1];
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
  Difference<Value> sum{};
  for (int j = 0; j <= n - 2; ++j) {
    numerator += numerators[j];
    Difference<Value> difference = values[j] - values[j + 1];
    sum += numerator * difference;
  }
  return sum / (denominator * step);
} else if (IMPL == 4) {
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
} else if (IMPL == 5) {
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
    sum += DoublePrecision<Difference<Value>>(numerator * difference);
  }
  return sum.value / (denominator * step);
} else { CHECK_EQ(IMPL, 6);
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
    sum += TwoProduct(numerator, difference);
  }
  return sum.value / (denominator * step);
}
}

}  // namespace internal_finite_difference

using internal_finite_difference::FiniteDifference;

}  // namespace numerics
}  // namespace principia
