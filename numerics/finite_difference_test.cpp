#include "numerics/finite_difference.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {
namespace numerics {

using ::testing::Each;
using ::testing::Lt;
using ::testing::Types;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_finite_difference;
using namespace principia::quantities::_quantities;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_statistics;

template<typename T>
class FiniteDifferenceTest : public ::testing::Test {
 protected:
  static constexpr int n = T::value;
};

template<int... values>
using IntegralConstants = Types<std::integral_constant<int, values>...>;
using UpTo9 = IntegralConstants<1, 2, 3, 4, 5, 6, 7, 8, 9>;

TYPED_TEST_SUITE(FiniteDifferenceTest, UpTo9);

// Polynomials of degree up to n - 1, differentiated exactly.
TYPED_TEST(FiniteDifferenceTest, LowDegreePolynomials) {
  constexpr int n = TestFixture::n;
  std::array<double, n> values;
  constexpr double h = 3;
  static constexpr std::int64_t max_bits =
      (std::array<std::int64_t, 9>{{0, 0, 5, 6, 7, 11, 13, 18, 21}})[n - 1];
  for (int p = 0; p < n; ++p) {
    for (int i = 0; i < n; ++i) {
      values[i] = 1729 + std::pow(i * h + π, p);
    }
    for (int j = 0; j < n; ++j) {
      EXPECT_THAT(
          FiniteDifference(values, /*step=*/h, /*offset=*/j),
          AlmostEquals(p * std::pow(j * h + π, p - 1),
                       0,
                       std::exp2(max_bits) - 1))
          << "with n = " << n << ", p = " << p << ", j = " << j;
    }
  }
}

// Polynomial of degree n, differentiated approximately, with convergence order
// n - 1.
TYPED_TEST(FiniteDifferenceTest, HighDegreePolynomial) {
  constexpr int n = TestFixture::n;
  std::array<double, n> values;
  for (int j = 0; j < n; ++j) {
    std::vector<double> log_steps;
    std::vector<double> log_errors;
    double const actual_derivative = n * std::pow(π, n - 1);
    for (double h = 1; h > 0x1p-3; h /= 2) {
      for (int i = 0; i < n; ++i) {
        values[i] = std::pow((i - j) * h + π, n);
      }
      log_steps.push_back(std::log2(h));
      log_errors.push_back(std::log2(
          RelativeError(FiniteDifference(values, /*step=*/h, /*offset=*/j),
                        actual_derivative)));
    }
    if constexpr (n == 1) {
      EXPECT_THAT(log_errors, Each(Infinity<double>));
    } else {
      EXPECT_THAT(AbsoluteError(n - 1, Slope(log_steps, log_errors)),
                  Lt(0.47))
          << "with n = " << n << ", j = " << j;
    }
  }
}

}  // namespace numerics
}  // namespace principia
