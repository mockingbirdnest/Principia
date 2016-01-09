#include "numerics/root_finders.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using quantities::Sqrt;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;

namespace numerics {

class RootFindersTest : public ::testing::Test {};

TEST_F(RootFindersTest, SquareRoots) {
  // Solving x * x == n for integers n.
  double const x_max = 10;
  double const n_max = x_max * x_max;
  for (double n = 1; n < n_max; ++n) {
    int evaluations = 0;
    auto const equation = [&n, &evaluations](double const x) {
      ++evaluations;
      return x * x - n;
    };
    EXPECT_THAT(Bisect(equation, 0.0, x_max), AlmostEquals(Sqrt(n), 0, 1));
    if (n == 25) {
      EXPECT_EQ(3, evaluations);
    } else {
      EXPECT_THAT(evaluations, AllOf(Ge(49), Le(58)));
    }
  }
}

}  // namespace numerics
}  // namespace principia
