
#include "testing_utilities/is_near.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/bipm.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace testing_utilities {

using quantities::Speed;
using quantities::bipm::Knot;
using quantities::si::Metre;
using quantities::si::Second;

class IsNearTest : public testing::Test {};

TEST_F(IsNearTest, Dimensionless) {
  double const y = e;
  EXPECT_THAT(y, IsNear(e, 1.0));
  EXPECT_THAT(y, Not(IsNear(3.0, 1.01)));
}

TEST_F(IsNearTest, Quantity) {
  Speed v = 1 * Knot;
  EXPECT_THAT(v, IsNear(0.514 * Metre / Second, 1.002));
  EXPECT_THAT(v, Not(IsNear(3.0 * Metre / Second, 1.01)));
}

TEST_F(IsNearTest, Negatives) {
  EXPECT_THAT(π - std::exp(π), IsNear(-20, 1.00001));
}

}  // namespace testing_utilities
}  // namespace principia
