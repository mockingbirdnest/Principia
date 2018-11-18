
#include "numerics/max_abs_normalized_associated_legendre_function.mathematica.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/combinatorics.hpp"
#include "numerics/legendre.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using quantities::Sqrt;
using testing_utilities::AlmostEquals;

class MaxAbsNormalizedAssociatedLegendreFunctionTest : public testing::Test {
 protected:
  static constexpr double ApproximateFactorial(std::int64_t const n) {
    double result = 1;
    for (std::int64_t i = 1; i <= n; ++i) {
      result *= i;
    }
    return result;
  }

  static constexpr double ApproximateDoubleFactorial(std::int64_t const n) {
    double result = 1;
    for (std::int64_t k = n; k > 0; k -= 2) {
      result *= k;
    }
    return result;
  }
};

TEST_F(MaxAbsNormalizedAssociatedLegendreFunctionTest, Polynomials) {
  for (int n = 0; n < MaxAbsNormalizedAssociatedLegendreFunction.rows; ++n) {
    EXPECT_THAT(MaxAbsNormalizedAssociatedLegendreFunction[n][0],
                AlmostEquals(LegendreNormalizationFactor[n][0], 0))
        << "n = " << n;
  }
}

TEST_F(MaxAbsNormalizedAssociatedLegendreFunctionTest, HighestOrder) {
  for (int n = 0; n < MaxAbsNormalizedAssociatedLegendreFunction.rows; ++n) {
    // TODO(phl): This should be
    //  LegendreNormalizationFactor[n][n] * DoubleFactorial(2 * n - 1),
    // but that NaNs.
    EXPECT_THAT(MaxAbsNormalizedAssociatedLegendreFunction[n][n],
                AlmostEquals(Sqrt((2 * n + 1) * (2 - (n == 0 ? 1 : 0)) /
                                  ApproximateFactorial(2 * n)) *
                                 ApproximateDoubleFactorial(2 * n - 1), 0, 1))
        << "n = " << n;
  }
}

}  // namespace numerics
}  // namespace principia
