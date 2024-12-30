#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"  // 🧙 For π.
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace functions {

using ::testing::AnyOf;
using namespace boost::multiprecision;
using namespace principia::functions::_multiprecision;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;

class StdAccuracyTest : public ::testing::Test {
 protected:
  static double ULPDistance(cpp_bin_float_50 const& actual,
                            cpp_bin_float_50 const& expected) {
    std::int64_t actual_exponent;
    std::int64_t expected_exponent;
    auto actual_mantissa = frexp(actual, &actual_exponent);
    auto expected_mantissa = frexp(expected, &expected_exponent);
    if (actual_exponent == expected_exponent) {
    } else if (actual_exponent == expected_exponent + 1) {
      --actual_exponent;
      actual_mantissa *= 2.0;
    } else if (actual_exponent == expected_exponent - 1) {
      ++actual_exponent;
      actual_mantissa /= 2.0;
    } else {
      LOG(FATAL) << actual_exponent << " " << actual_mantissa << " "
                 << expected_exponent << " " << expected_mantissa;
    }
    return std::abs(static_cast<double>(actual_mantissa - expected_mantissa)) *
           (1LL << std::numeric_limits<double>::digits);
  }
};

#if !_DEBUG

TEST_F(StdAccuracyTest, SinCos) {
  // Random values without argument reduction.
  {
    std::mt19937_64 random(42);
    std::uniform_real_distribution<> angle_distribution(0, π / 4);
    double max_sin_ulp_distance = 0;
    double max_cos_ulp_distance = 0;
    for (int i = 0; i < 1e5; ++i) {
      double const α = angle_distribution(random);
      max_sin_ulp_distance =
          std::max(max_sin_ulp_distance, ULPDistance(std::sin(α), Sin(α)));
      max_cos_ulp_distance =
          std::max(max_cos_ulp_distance, ULPDistance(std::cos(α), Cos(α)));
    }
    EXPECT_THAT(max_sin_ulp_distance,
                AnyOf(IsNear(0.727_(1)),    // Windows.
                      IsNear(0.654_(1)),    // macOS.
                      IsNear(0.513_(1))));  // Ubuntu.
    EXPECT_THAT(max_cos_ulp_distance,
                AnyOf(IsNear(0.834_(1)),    // Windows, macOS.
                      IsNear(0.502_(1))));  // Ubuntu.
  }

  // Hardest argument reduction, [Mul+10] table 11.1.
  {
    double const x = 0x16ac5b262ca1ffp797;
    EXPECT_THAT(ULPDistance(std::sin(x), Sin(x)), IsNear(9.89e-22_(1)));
    EXPECT_THAT(ULPDistance(std::cos(x), Cos(x)),
                AnyOf(IsNear(0.0454_(1)),  // Windows, macOS.
                      IsNear(7.95_(1))));  // Linux.
  }
}

#endif

}  // namespace functions
}  // namespace principia
