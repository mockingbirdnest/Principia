#include <algorithm>
#include <limits>
#include <random>

#include "core-math/cos.h"
#include "core-math/sin.h"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"  // 🧙 For π.
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace functions {

using namespace boost::multiprecision;
using namespace principia::functions::_multiprecision;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;

class CoreMathAccuracyTest : public ::testing::Test {
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

TEST_F(CoreMathAccuracyTest, SinCos) {
  // Random values without argument reduction.
  {
    std::mt19937_64 random(42);
    std::uniform_real_distribution<> angle_distribution(0, π / 4);
    double max_sin_ulp_distance = 0;
    double max_cos_ulp_distance = 0;
    for (int i = 0; i < 1e5; ++i) {
      double const α = angle_distribution(random);
      max_sin_ulp_distance =
          std::max(max_sin_ulp_distance, ULPDistance(cr_sin(α), Sin(α)));
      max_cos_ulp_distance =
          std::max(max_cos_ulp_distance, ULPDistance(cr_cos(α), Cos(α)));
    }
    EXPECT_THAT(max_sin_ulp_distance, IsNear(0.49999_(1)));
    EXPECT_THAT(max_cos_ulp_distance, IsNear(0.49999_(1)));
  }

  // Hardest argument reduction, [Mul+10] table 11.1.
  {
    double const x = 0x16ac5b262ca1ffp797;
    EXPECT_THAT(ULPDistance(cr_sin(x), Sin(x)), IsNear(9.89e-22_(1)));
    EXPECT_THAT(ULPDistance(cr_cos(x), Cos(x)), IsNear(0.0454_(1)));
  }
}

#endif

}  // namespace functions
}  // namespace principia
