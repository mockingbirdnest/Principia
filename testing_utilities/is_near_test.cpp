
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
using quantities::bipm::NauticalMile;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Second;

class IsNearTest : public testing::Test {};

TEST_F(IsNearTest, Dimensionless) {
  double const y = e;
  EXPECT_THAT(y, IsNear(2.718_⑴));
  EXPECT_THAT(y, Not(IsNear(3.0_⑵)));
}

TEST_F(IsNearTest, Quantity) {
  Speed v = 1 * Knot;
  EXPECT_THAT(v, IsNear(0.514_⑴ * Metre / Second));
  EXPECT_THAT(v, Not(IsNear(3.0_⑼ * Metre / Second)));
  EXPECT_THAT(v, IsNear(1.0_⑴ * NauticalMile / Hour));
  EXPECT_THAT(v, Not(IsNear(1.2_⑴ * NauticalMile / Hour)));
}

TEST_F(IsNearTest, Negatives) {
  EXPECT_THAT(π - std::exp(π), IsNear(-20.000_⑴));
}

}  // namespace testing_utilities
}  // namespace principia
