#include <cfloat>
#include <cmath>

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace testing_utilities {

using si::Metre;
using testing::Eq;
using testing::Ne;

namespace {
struct World;
}  // namespace

class NumericsTest : public testing::Test {};

// The smallest positive double, a denormal.
double const SmallestPositive = DBL_MIN * DBL_EPSILON;

TEST_F(NumericsTest, ULPs) {
  EXPECT_THAT(ULPDistance(1, 1), Eq(0));
  EXPECT_THAT(ULPDistance(+0.0, +0.0), Eq(0));
  EXPECT_THAT(ULPDistance(+0.0, -0.0), Eq(0));
  // DBL_MIN is the smallest positive normalised number.
  // 52 bits of mantissa stand between it and 0, in the form of denormals.
  EXPECT_THAT(ULPDistance(+0.0, DBL_MIN), Eq(std::pow(2, DBL_MANT_DIG - 1)));
  EXPECT_THAT(ULPDistance(-0.0, DBL_MIN), Eq(std::pow(2, DBL_MANT_DIG - 1)));
  EXPECT_THAT(ULPDistance(+0.0, -SmallestPositive), Eq(1));
  EXPECT_THAT(ULPDistance(-0.0, -SmallestPositive), Eq(1));
  EXPECT_THAT(ULPDistance(-1, 1), Ne(0));
  EXPECT_THAT(ULPDistance(-1, 1), 2 * ULPDistance(0, 1));
}

TEST_F(NumericsTest, AbsoluteError) {
  auto const double_abs = [](double const x) { return std::abs(x); };
  EXPECT_THAT(AbsoluteError(1., 1., double_abs), Eq(0.));
  EXPECT_THAT(AbsoluteError(1., 2., double_abs), Eq(1.));
  EXPECT_THAT(AbsoluteError(1., 0., double_abs), Eq(1.));
  EXPECT_THAT(AbsoluteError(1, 1), Eq(0));
  EXPECT_THAT(AbsoluteError(1, 2), Eq(1));
  EXPECT_THAT(AbsoluteError(1, 0), Eq(1));
  EXPECT_THAT(AbsoluteError(1 * Metre, 1 * Metre), Eq(0 * Metre));
  EXPECT_THAT(AbsoluteError(1 * Metre, 2 * Metre), Eq(1 * Metre));
  EXPECT_THAT(AbsoluteError(1 * Metre, 0 * Metre), Eq(1 * Metre));
}

}  // namespace testing_utilities
}  // namespace principia
