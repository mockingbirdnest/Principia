
#include "numerics/fast_sin_cos_2π.hpp"

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

class FastSinCos2πTest : public ::testing::Test {
 protected:
#if defined(_DEBUG)
  static constexpr int iterations_ = 1e7;
#else
  static constexpr int iterations_ = 1e8;
#endif
};

// Check that we obtain the right result for special angles.
// TODO(phl): We don't.  Improve this.
TEST_F(FastSinCos2πTest, SpecialValues) {
  double sin;
  double cos;
  FastSinCos2π(0.0, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(0.0, 0));
  EXPECT_THAT(cos, AlmostEquals(1.0, 0));
  FastSinCos2π(0.25, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(1.0, 0));
  EXPECT_THAT(cos, AlmostEquals(0.0, 0));
  FastSinCos2π(0.5, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(0.0, 0));
  EXPECT_THAT(cos, AlmostEquals(-1.0, 0));
  FastSinCos2π(0.75, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(-1.0, 0));
  EXPECT_THAT(cos, AlmostEquals(0.0, 0));
  FastSinCos2π(1.0, sin, cos);
  EXPECT_THAT(sin, AlmostEquals(0.0, 0));
  EXPECT_THAT(cos, AlmostEquals(1.0, 0));
}

// Check that, for arguments in [-1, 1], the errors are of the expected size and
// at the expected location.
TEST_F(FastSinCos2πTest, Random1) {
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
    FastSinCos2π(x, sin, cos);
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

// Arguments in the range [-1e6, 1e6 + 1] to test argument reduction.
TEST_F(FastSinCos2πTest, Random1000) {
  std::mt19937_64 random(42);
  std::uniform_int_distribution<> const integer_distribution(-1e6, 1e6);
  std::uniform_real_distribution<> const fractional_distribution(0.0, 1.0);
  double max_sin_error = 0.0;
  double max_cos_error = 0.0;
  double max_sin_error_x = 0.0;
  double max_cos_error_x = 0.0;
  for (int i = 0; i < iterations_; ++i) {
    double const integer = integer_distribution(random);
    double const fractional = fractional_distribution(random);
    double const x = integer + fractional;
    double const x_reduced = x - integer;
    double sin;
    double cos;
    FastSinCos2π(x, sin, cos);
    double sin_reduced;
    double cos_reduced;
    FastSinCos2π(x_reduced, sin_reduced, cos_reduced);

    // Smart expectations are too slow, so let's be dumb.
    EXPECT_TRUE(sin_reduced == sin);
    EXPECT_TRUE(cos_reduced == cos);
  }
}

// Check that the functions are reasonably monotonic.
TEST_F(FastSinCos2πTest, Monotonicity) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> const distribution(0.0, 0.25);
  double max_sin_error = 0.0;
  double max_cos_error = 0.0;
  double max_sin_error_x = 0.0;
  double max_cos_error_x = 0.0;
  for (int i = 0; i < iterations_; ++i) {
    double const x = distribution(random);
    double const next_x = std::nexttoward(x, 1.0);
    double sin_2π_x;
    double cos_2π_x;
    FastSinCos2π(x, sin_2π_x, cos_2π_x);
    double sin_2π_next_x;
    double cos_2π_next_x;
    FastSinCos2π(next_x, sin_2π_next_x, cos_2π_next_x);

    if (sin_2π_next_x < sin_2π_x) {
      max_sin_error = std::max(max_sin_error, sin_2π_x - sin_2π_next_x);
      max_sin_error_x = x;
    }
    if (cos_2π_next_x > cos_2π_x) {
      max_cos_error = std::max(max_cos_error, cos_2π_next_x - cos_2π_x);
      max_cos_error_x = x;
    }
  }
  EXPECT_LT(max_sin_error, 4e-16);
  EXPECT_LT(max_cos_error, 4e-16);
}

}  // namespace numerics
}  // namespace principia
