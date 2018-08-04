
#include "numerics/fast_sin_cos_cycle.hpp"

#include <algorithm>
#include <random>

#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {

using testing_utilities::IsNear;

namespace numerics {

class FastSinCosCycleTest : public ::testing::Test {};

// Check that, for arguments in [-1, 1], the errors are of the expected size and
// at the expected location.
TEST_F(FastSinCosCycleTest, Random1) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> const distribution(-1.0, 1.0);
  double max_sin_error = 0.0;
  double max_cos_error = 0.0;
  double max_sin_error_x = 0.0;
  double max_cos_error_x = 0.0;
  for (int i = 0; i < 1e8; ++i) {
    double const x = distribution(random);
    double sin;
    double cos;
    FastSinCosCycle(x, sin, cos);
    double const sin_error = std::abs(std::sin(2 * π * x) - sin);
    if (sin_error > max_sin_error) {
      max_sin_error = sin_error;
      max_sin_error_x = x;
    }
    double const cos_error = std::abs(std::cos(2 * π * x) - cos);
    if (cos_error > max_cos_error) {
      max_cos_error = cos_error;
      max_cos_error_x = x;
    }
  }
  // These numbers come from the Mathematica minimax computation.
  EXPECT_LT(max_sin_error, 5.90e-7);
  EXPECT_LT(max_cos_error, 5.28e-8);
  EXPECT_THAT(max_sin_error_x, IsNear(-1.0 + 1.0 / 23.0));
  EXPECT_THAT(max_cos_error_x, IsNear(1.0 - 1.0 / 15.0));
}

}  // namespace numerics
}  // namespace principia
