#include "testing_utilities/optimization_test_functions.hpp"

#include "gtest/gtest.h"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace testing_utilities {

using ::testing::ElementsAre;
using ::testing::Eq;

TEST(OptimizationTestFunctionsTest, Branin) {
  {
    double const x₁ = -π;
    double const x₂ = 12.275;
    EXPECT_THAT(Branin(x₁, x₂), IsNear(0.397887_(1)));
    EXPECT_THAT(𝛁Branin(x₁, x₂),
                ElementsAre(IsNear(1.2e-15_(1)), Eq(0)));
  }
  {
    double const x₁ = π;
    double const x₂ = 2.275;
    EXPECT_THAT(Branin(x₁, x₂), IsNear(0.397887_(1)));
    EXPECT_THAT(𝛁Branin(x₁, x₂),
                ElementsAre(IsNear(-4.8e-16_(1)), IsNear(8.9e-16_(1))));
  }
  {
    double const x₁ = 9.42478;
    double const x₂ = 2.475;
    EXPECT_THAT(Branin(x₁, x₂), IsNear(0.397887_(1)));
    EXPECT_THAT(𝛁Branin(x₁, x₂),
                ElementsAre(IsNear(2.2e-5_(1)), IsNear(-3.4e-6_(1))));
  }
  {
    double const x₁ = 0.5;
    double const x₂ = -0.3;
    EXPECT_THAT(Branin(x₁, x₂), IsNear(49.0797_(1)));
    EXPECT_THAT(𝛁Branin(x₁, x₂),
                ElementsAre(IsNear(-20.7963_(1)), IsNear(-11.073_(1))));
  }
}

TEST(OptimizationTestFunctionsTest, GoldsteinPrice) {
  {
    double const x₁ = 0;
    double const x₂ = -1;
    EXPECT_THAT(GoldsteinPrice(x₁, x₂), Eq(3));
    EXPECT_THAT(𝛁GoldsteinPrice(x₁, x₂),
                ElementsAre(Eq(0), Eq(0)));
  }
  {
    double const x₁ = 0.5;
    double const x₂ = -0.3;
    EXPECT_THAT(GoldsteinPrice(x₁, x₂), IsNear(596.161_(1)));
    EXPECT_THAT(𝛁GoldsteinPrice(x₁, x₂),
                ElementsAre(IsNear(-601.51_(1)), IsNear(2163.65_(1))));
  }
}

TEST(OptimizationTestFunctionsTest, Hartmann3) {
  {
    double const x₁ = 0.114589;
    double const x₂ = 0.555649;
    double const x₃ = 0.852547;
    EXPECT_THAT(Hartmann3(x₁, x₂, x₃), IsNear(-3.86278_(1)));
    EXPECT_THAT(𝛁Hartmann3(x₁, x₂, x₃),
                ElementsAre(IsNear(5.06186e-8_(1)),
                            IsNear(-1.26603e-5_(1)),
                            IsNear(-7.74866e-05_(1))));
  }
  {
    double const x₁ = 0.5;
    double const x₂ = -0.3;
    double const x₃ = 0.8;
    EXPECT_THAT(Hartmann3(x₁, x₂, x₃), IsNear(-0.00585298_(1)));
    EXPECT_THAT(𝛁Hartmann3(x₁, x₂, x₃),
                ElementsAre(IsNear(0.000101777_(1)),
                            IsNear(-0.0894713_(1)),
                            IsNear(0.0113956_(1))));
  }
}

}  // namespace testing_utilities
}  // namespace principia
