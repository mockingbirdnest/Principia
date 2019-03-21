#pragma once

#include "numerics/finite_difference.hpp"

#include "numerics/double_precision.hpp"
#include "numerics/finite_difference.mathematica.h"

namespace principia {
namespace numerics {
namespace internal_finite_difference {

// Represents a linear combination of values with coefficients |numerators[i]|
// divided by the common |denominator|.
// The sum of the entries in |numerators| must be 0.

constexpr int max_order = 16;

constexpr FixedLowerTriangularMatrix<double, max_order + 1> forward_differences{
    {{0, -1, 1, -3, 4, -1, -11, 18, -9, 2, -25, 48, -36, 16, -3, -137, \
300, -300, 200, -75, 12, -147, 360, -450, 400, -225, 72, -10, -1089, \
2940, -4410, 4900, -3675, 1764, -490, 60, -2283, 6720, -11760, 15680, \
-14700, 9408, -3920, 960, -105, -7129, 22680, -45360, 70560, -79380, \
63504, -35280, 12960, -2835, 280, -7381, 25200, -56700, 100800, \
-132300, 127008, -88200, 43200, -14175, 2800, -252, -83711, 304920, \
-762300, 1524600, -2286900, 2561328, -2134440, 1306800, -571725, \
169400, -30492, 2520, -86021, 332640, -914760, 2032800, -3430350, \
4390848, -4268880, 3136320, -1715175, 677600, -182952, 30240, -2310, \
-1145993, 4684680, -14054040, 34354320, -64414350, 92756664, \
-103062960, 88339680, -57972915, 28628600, -10306296, 2555280, \
-390390, 27720, -1171733, 5045040, -16396380, 43723680, -90180090, \
144288144, -180360180, 176679360, -135270135, 80160080, -36072036, \
11924640, -2732730, 388080, -25740, -1195757, 5405400, -18918900, \
54654600, -122972850, 216432216, -300600300, 331273800, -289864575, \
200400200, -108216108, 44717400, -13663650, 2910600, -386100, 24024, \
-2436559, 11531520, -43243200, 134534400, -327927600, 629620992, \
-961920960, 1177862400, -1159458300, 916115200, -577152576, \
286191360, -109309200, 31046400, -6177600, 768768, -45045}}};

constexpr std::array<double, max_order + 1> denominators{{1,
                                                          1,
                                                          2,
                                                          6,
                                                          12,
                                                          60,
                                                          60,
                                                          420,
                                                          840,
                                                          2520,
                                                          2520,
                                                          27720,
                                                          27720,
                                                          360360,
                                                          360360,
                                                          360360,
                                                          720720}};

static int IMPL;

template<typename Value, typename Argument, int n>
Derivative<Value, Argument> FiniteDifference(
    FixedVector<Value, n> const& values,
    Argument const& step,
    int offset) {
  CHECK_EQ(offset, 0);
  constexpr double const* numerators = forward_differences[n - 1];
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
