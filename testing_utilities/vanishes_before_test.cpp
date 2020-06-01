
#include "testing_utilities/vanishes_before.hpp"

#include <limits>
#include <sstream>

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/bipm.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using quantities::Speed;
using quantities::bipm::Knot;
using ::testing::Ne;
namespace si = quantities::si;

namespace testing_utilities {

class VanishesBeforeTest : public testing::Test {};

TEST_F(VanishesBeforeTest, Dimensionless) {
  double const y = 3000.0 * std::numeric_limits<double>::epsilon();
  EXPECT_THAT(y, VanishesBefore(1000.0, 6));
  EXPECT_THAT(2 * y, Not(VanishesBefore(1000.0, 6)));
  double const δy = e / 100.0;
  double e_accumulated = 0.0;
  for (int i = 1; i <= 100.0; ++i) {
    e_accumulated += δy;
  }
  EXPECT_THAT(e_accumulated, Ne(e));
  EXPECT_THAT(e_accumulated - e, Not(VanishesBefore(e, 4)));
  EXPECT_THAT(e_accumulated - e, VanishesBefore(e, 1));
}

TEST_F(VanishesBeforeTest, Quantity) {
  Speed v1 = 1 * Knot;
  Speed const v2 = 3 * v1 * std::numeric_limits<double>::epsilon();
  EXPECT_THAT(v2, VanishesBefore(v1, 3));
  EXPECT_THAT(2 * v2, Not(VanishesBefore(v1, 3)));
  Speed const δv = v1 / 100.0;
  Speed v_accumulated;
  // If the upper bound is an int, the compiler adds 10 times 10 * δv and this
  // results in no error in the final result.
  for (int i = 1; i <= 100.0; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated - v1, Not(VanishesBefore(v1, 8)));
  EXPECT_THAT(v_accumulated - v1, VanishesBefore(v1, 4));
}

TEST_F(VanishesBeforeTest, Describe) {
  Speed v1 = 1 * si::Unit<Speed>;
  {
    std::ostringstream out;
    VanishesBefore(v1, 2, 6).impl().DescribeTo(&out);
    EXPECT_EQ("vanishes before +1.00000000000000000e+00 m s^-1 "
              "to within 2 to 6 ULPs",
              out.str());
  }
  {
    std::ostringstream out;
    VanishesBefore(v1, 2, 6).impl().DescribeNegationTo(&out);
    EXPECT_EQ("does not vanish before +1.00000000000000000e+00 m s^-1 "
              "to within 2 to 6 ULP",
              out.str());
  }
}

}  // namespace testing_utilities
}  // namespace principia
