#include <cfloat>
#include <cmath>

#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace testing_utilities {

using quantities::Dimensionless;
using quantities::Sqrt;
using geometry::R3Element;
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

TEST_F(NumericsTest, DoubleAbsoluteError) {
  auto const double_abs = [](double const x) { return std::abs(x); };
  EXPECT_THAT(AbsoluteError(1., 1., double_abs), Eq(0.));
  EXPECT_THAT(AbsoluteError(1., 2., double_abs), Eq(1.));
  EXPECT_THAT(AbsoluteError(1., 0., double_abs), Eq(1.));
}

TEST_F(NumericsTest, DimensionlessAbsoluteError) {
  EXPECT_THAT(AbsoluteError(1, 1), Eq(0));
  EXPECT_THAT(AbsoluteError(1, 2), Eq(1));
  EXPECT_THAT(AbsoluteError(1, 0), Eq(1));
}

TEST_F(NumericsTest, DimensionfulAbsoluteError) {
  EXPECT_THAT(AbsoluteError(1 * Metre, 1 * Metre), Eq(0 * Metre));
  EXPECT_THAT(AbsoluteError(1 * Metre, 2 * Metre), Eq(1 * Metre));
  EXPECT_THAT(AbsoluteError(1 * Metre, 0 * Metre), Eq(1 * Metre));
}

TEST_F(NumericsTest, R3ElementAbsoluteError) {
  R3Element<Dimensionless> const i = {1, 0, 0};
  R3Element<Dimensionless> const j = {0, 1, 0};
  R3Element<Dimensionless> const k = {0, 0, 1};
  EXPECT_THAT(AbsoluteError(i + j + k, i + j + k), Eq(0));
  EXPECT_THAT(AbsoluteError(i + j, i + j + k), Eq(1));
  EXPECT_THAT(AbsoluteError(i, i + j + k), Eq(Sqrt(2)));
}

}  // namespace testing_utilities
}  // namespace principia
