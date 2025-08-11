
#include "numerics/max_abs_normalized_associated_legendre_function.mathematica.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/combinatorics.hpp"
#include "numerics/legendre_normalization_factor.mathematica.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using ::testing::Lt;
using namespace principia::numerics::_legendre_normalization_factor;
using namespace principia::numerics::_max_abs_normalized_associated_legendre_function;  // NOLINT
using namespace principia::numerics::_elementary_functions;
using namespace principia::testing_utilities::_almost_equals;

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
  for (int n = 0; n < MaxAbsNormalizedAssociatedLegendreFunction.rows(); ++n) {
    EXPECT_THAT(MaxAbsNormalizedAssociatedLegendreFunction(n, 0),
                AlmostEquals(LegendreNormalizationFactor(n, 0), 0))
        << "n = " << n;
  }
}

TEST_F(MaxAbsNormalizedAssociatedLegendreFunctionTest, HighestOrder) {
  for (int n = 0; n < MaxAbsNormalizedAssociatedLegendreFunction.rows(); ++n) {
    EXPECT_THAT(MaxAbsNormalizedAssociatedLegendreFunction(n, n),
                AlmostEquals(LegendreNormalizationFactor(n, n) *
                                 ApproximateDoubleFactorial(2 * n - 1),
                             0, 3))
        << "n = " << n;
  }
}

}  // namespace numerics
}  // namespace principia
