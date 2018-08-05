
#include "numerics/fast_sin_cos_cycle.hpp"

#include <algorithm>
#include <random>

#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::VanishesBefore;

namespace numerics {

class FastSinCosCycleTest : public ::testing::Test {
 protected:
#if defined(_DEBUG)
  static constexpr int iterations_ = 1e7;
#else
  static constexpr int iterations_ = 1e8;
#endif
};

// Check that we obtain the right result for special angles.
// TODO(phl): We don't.  Improve this.
TEST_F(FastSinCosCycleTest, SpecialValues) {
  double sin;
  double cos;
  FastSinCosCycle(0.0, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(0.0, 0));
  EXPECT_THAT(cos, AlmostEquals(1.0, 0));
  FastSinCosCycle(0.25, sin, cos);
  EXPECT_THAT(sin, IsNear(1.0, 1.00001));
  EXPECT_THAT(cos, VanishesBefore(1.0, 0, 3e8));
  FastSinCosCycle(0.5, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(0.0, 0));
  EXPECT_THAT(cos, AlmostEquals(-1.0, 0));
  FastSinCosCycle(0.75, sin, cos);
  EXPECT_THAT(sin, IsNear(-1.0));
  EXPECT_THAT(cos, VanishesBefore(1.0, 0, 3e8));
  FastSinCosCycle(1.0, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(0.0, 0));
  EXPECT_THAT(cos, AlmostEquals(1.0, 0));
}

// Check that, for arguments in [-1, 1], the errors are of the expected size and
// at the expected location.
TEST_F(FastSinCosCycleTest, Random1) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> const distribution(-1.0, 1.0);
  double max_sin_error = 0.0;
  double max_cos_error = 0.0;
  double max_sin_error_x = 0.0;
  double max_cos_error_x = 0.0;
  for (int i = 0; i < iterations_; ++i) {
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

// Same as above, but in the range [-1e6, 1e6] to test argument reduction.
TEST_F(FastSinCosCycleTest, Random1000) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> const distribution(-1.0e6, 1.0e6);
  double max_sin_error = 0.0;
  double max_cos_error = 0.0;
  double max_sin_error_x = 0.0;
  double max_cos_error_x = 0.0;
  for (int i = 0; i < iterations_; ++i) {
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
  // Slightly larger errors here because of argument reduction.
  EXPECT_LT(max_sin_error, 5.90e-7);
  EXPECT_LT(max_cos_error, 5.35e-8);
}

// Check that the functions are reasonably monotonic.
TEST_F(FastSinCosCycleTest, Monotonicity) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> const distribution(0.0, 0.25);
  double max_sin_error = 0.0;
  double max_cos_error = 0.0;
  double max_sin_error_x = 0.0;
  double max_cos_error_x = 0.0;
  bool found_anomaly = false;
  for (int i = 0; i < iterations_; ++i) {
    double const x = distribution(random);
    double const next_x = std::nexttoward(x, 1.0);
    double sin_x;
    double cos_x;
    FastSinCosCycle(x, sin_x, cos_x);
    double sin_next_x;
    double cos_next_x;
    FastSinCosCycle(x, sin_next_x, cos_next_x);

    if (sin_next_x < sin_x) {
      found_anomaly = true;
      max_sin_error = std::max(max_sin_error, sin_x - sin_next_x);
      max_sin_error_x = x;
    }
    if (cos_next_x > cos_x) {
      found_anomaly = true;
      max_cos_error = std::max(max_cos_error, cos_next_x - cos_x);
      max_cos_error_x = x;
    }
  }
  EXPECT_FALSE(found_anomaly)
      << "sin error: " << max_sin_error << " at " << max_sin_error_x
      << ", cos error: " << max_cos_error << " at " << max_cos_error_x;
}

}  // namespace numerics
}  // namespace principia
