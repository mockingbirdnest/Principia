#include "testing_utilities/optimization_test_functions.hpp"

#include "gtest/gtest.h"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace testing_utilities {

TEST(OptimizationTestFunctionsTest, GoldsteinPrice) {
  double const x₁ = 0.5;
  double const x₂ = -0.3;
  EXPECT_THAT(GoldsteinPrice(x₁, x₂), IsNear(596.161_(1)));
  EXPECT_THAT(GradGoldsteinPrice(x₁, x₂),
              ElementsAre(IsNear(-601.51_(1)), IsNear(2163.65_(1))));
}

}  // namespace testing_utilities
}  // namespace principia
