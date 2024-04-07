#include <cmath>
#include <random>

#include "functions/multiprecision.hpp"
#include "gtest/gtest.h"
#include "numerics/ulp_distance.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {

using namespace boost::multiprecision;
using namespace principia::numerics::_ulp_distance;
using namespace principia::quantities::_si;

class StdAccuracyTest : public ::testing::Test {};

#if !_DEBUG

TEST_F(StdAccuracyTest, SinCos) {
  // Random values without argument reduction.
  {
    std::mt19937_64 random(42);
    std::uniform_real_distribution<> angle_distribution(0, π / 4);
    std::int64_t max_sin_ulp_distance = 0;
    std::int64_t max_cos_ulp_distance = 0;
    for (int i = 0; i < 1e5; ++i) {
      double const α = angle_distribution(random);
      max_sin_ulp_distance =
          std::max(max_sin_ulp_distance,
                   ULPDistance(std::sin(α), static_cast<double>(Sin(α))));
      max_cos_ulp_distance =
          std::max(max_cos_ulp_distance,
                   ULPDistance(std::cos(α), static_cast<double>(Cos(α))));
    }
    EXPECT_EQ(1, max_sin_ulp_distance);
    EXPECT_EQ(1, max_cos_ulp_distance);
  }

  // Hardest argument reduction, [Mul+10] table 11.1.
  {
    double const x = 0x16ac5b262ca1ffp797;
    EXPECT_EQ(
        0,
        ULPDistance(std::sin(x), static_cast<double>(Sin(x))));
    EXPECT_EQ(
        0,
        ULPDistance(std::cos(x), static_cast<double>(Cos(x))));
  }
}

#endif

}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
