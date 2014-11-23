#include "testing_utilities/numerics.hpp"

#include <cfloat>
#include <cmath>
#include <limits>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

using quantities::Sqrt;
using geometry::Bivector;
using geometry::R3Element;
using geometry::Trivector;
using geometry::Vector;
using si::Metre;
using testing::AllOf;
using testing::Eq;
using testing::Gt;
using testing::Lt;
using testing::Ne;

namespace {
struct World;
}  // namespace

class NumericsTest : public testing::Test {
 protected:
  struct World;

  R3Element<double> const i_ = {1, 0, 0};
  R3Element<double> const j_ = {0, 1, 0};
  R3Element<double> const k_ = {0, 0, 1};
};

double DoubleAbs(const double x) {
  return std::abs(x);
}

// The smallest positive double, a denormal.
double const SmallestPositive =
    DBL_MIN * std::numeric_limits<double>::epsilon();

TEST_F(NumericsTest, ULPs) {
  EXPECT_THAT(ULPDistance(1, 1), Eq(0));
  EXPECT_THAT(ULPDistance(+0.0, +0.0), Eq(0));
  EXPECT_THAT(ULPDistance(+0.0, -0.0), Eq(0));
  // DBL_MIN is the smallest positive normalized number.
  // 52 bits of mantissa stand between it and 0, in the form of denormals.
  EXPECT_THAT(ULPDistance(+0.0, DBL_MIN), Eq(std::pow(2, DBL_MANT_DIG - 1)));
  EXPECT_THAT(ULPDistance(-0.0, DBL_MIN), Eq(std::pow(2, DBL_MANT_DIG - 1)));
  EXPECT_THAT(ULPDistance(+0.0, -SmallestPositive), Eq(1));
  EXPECT_THAT(ULPDistance(-0.0, -SmallestPositive), Eq(1));
  EXPECT_THAT(ULPDistance(-1, 1), Ne(0));
  EXPECT_THAT(ULPDistance(-1, 1), Eq(2 * ULPDistance(0, 1)));
}

TEST_F(NumericsTest, DoubleAbsoluteError) {
  EXPECT_THAT(AbsoluteError(1., 1., DoubleAbs), Eq(0.));
  EXPECT_THAT(AbsoluteError(1., 2., DoubleAbs), Eq(1.));
  EXPECT_THAT(AbsoluteError(1., 0., DoubleAbs), Eq(1.));
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
  EXPECT_THAT(AbsoluteError(i_ + j_ + k_, i_ + j_ + k_), Eq(0));
  EXPECT_THAT(AbsoluteError(i_ + j_, i_ + j_ + k_), Eq(1));
  EXPECT_THAT(AbsoluteError(i_, i_ + j_ + k_), Eq(Sqrt(2)));
}

TEST_F(NumericsTest, VectorAbsoluteError) {
  EXPECT_THAT(
      AbsoluteError(Vector<double, World>(i_ + j_ + k_),
                    Vector<double, World>(i_ + j_ + k_)),
      Eq(0));
  EXPECT_THAT(
      AbsoluteError(Vector<double, World>(i_ + j_),
                    Vector<double, World>(i_ + j_ + k_)),
      Eq(1));
  EXPECT_THAT(
      AbsoluteError(Vector<double, World>(i_),
                    Vector<double, World>(i_ + j_ + k_)),
      Eq(Sqrt(2)));
}

TEST_F(NumericsTest, BivectorAbsoluteError) {
  EXPECT_THAT(
      AbsoluteError(Bivector<double, World>(i_ + j_ + k_),
                    Bivector<double, World>(i_ + j_ + k_)),
      Eq(0));
  EXPECT_THAT(
      AbsoluteError(Bivector<double, World>(i_ + j_),
                    Bivector<double, World>(i_ + j_ + k_)),
      Eq(1));
  EXPECT_THAT(
      AbsoluteError(Bivector<double, World>(i_),
                    Bivector<double, World>(i_ + j_ + k_)),
      Eq(Sqrt(2)));
}

TEST_F(NumericsTest, TrivectorAbsoluteError) {
  EXPECT_THAT(
      AbsoluteError(Trivector<double, World>(1),
                    Trivector<double, World>(1)),
      Eq(0));
  EXPECT_THAT(
      AbsoluteError(Trivector<double, World>(1.0),
                    Trivector<double, World>(2)),
      Eq(1));
  EXPECT_THAT(
      AbsoluteError(Trivector<double, World>(1),
                    Trivector<double, World>(0)),
      Eq(1));
}

TEST_F(NumericsTest, DoubleRelativeError) {
  EXPECT_THAT(RelativeError(42.0, 42.0, DoubleAbs), Eq(0));
  EXPECT_THAT(RelativeError(1.0, -1.0, DoubleAbs), Eq(2));
  EXPECT_THAT(RelativeError(2.0, 1.0, DoubleAbs), Eq(0.5));
  EXPECT_THAT(RelativeError(1.0, 2.0, DoubleAbs), Eq(1));
  EXPECT_THAT(RelativeError(42.0, 6.0 * 9.0, DoubleAbs),
              AllOf(Gt(0.28), Lt(0.29)));
}

TEST_F(NumericsTest, DimensionlessRelativeError) {
  EXPECT_THAT(RelativeError(42.0, 42.0), Eq(0));
  EXPECT_THAT(RelativeError(1.0, -1.0), Eq(2));
  EXPECT_THAT(RelativeError(2.0, 1.0), Eq(0.5));
  EXPECT_THAT(RelativeError(1.0, 2.0), Eq(1));
  EXPECT_THAT(RelativeError(42.0, 6.0 * 9.0), AllOf(Gt(0.28), Lt(0.29)));
}

TEST_F(NumericsTest, DimensionfulRelativeError) {
  EXPECT_THAT(RelativeError(42 * Metre, 42 * Metre), Eq(0));
  EXPECT_THAT(RelativeError(1 * Metre, -1 * Metre), Eq(2));
  EXPECT_THAT(RelativeError(2 * Metre, 1 * Metre), Eq(0.5));
  EXPECT_THAT(RelativeError(1 * Metre, 2 * Metre), Eq(1));
  EXPECT_THAT(RelativeError(42 * Metre, 6 * 9 * Metre),
              AllOf(Gt(0.28), Lt(0.29)));
}

TEST_F(NumericsTest, R3ElementRelativeError) {
  EXPECT_THAT(RelativeError(i_ + j_ + k_, i_ + j_ + k_), Eq(0));
  EXPECT_THAT(RelativeError(i_ + j_, i_ + j_ + k_), AlmostEquals(Sqrt(0.5)));
  EXPECT_THAT(RelativeError(i_, i_ + j_ + k_), Eq(Sqrt(2)));
}

TEST_F(NumericsTest, VectorRelativeError) {
  EXPECT_THAT(
      RelativeError(Vector<double, World>(i_ + j_ + k_),
                    Vector<double, World>(i_ + j_ + k_)),
      Eq(0));
  EXPECT_THAT(
      RelativeError(Vector<double, World>(i_ + j_),
                    Vector<double, World>(i_ + j_ + k_)),
      AlmostEquals(Sqrt(0.5)));
  EXPECT_THAT(
      RelativeError(Vector<double, World>(i_),
                    Vector<double, World>(i_ + j_ + k_)),
      Eq(Sqrt(2)));
}

TEST_F(NumericsTest, BivectorRelativeError) {
  EXPECT_THAT(
      RelativeError(Bivector<double, World>(i_ + j_ + k_),
                    Bivector<double, World>(i_ + j_ + k_)),
      Eq(0));
  EXPECT_THAT(
      RelativeError(Bivector<double, World>(i_ + j_),
                    Bivector<double, World>(i_ + j_ + k_)),
      AlmostEquals(Sqrt(0.5)));
  EXPECT_THAT(
      RelativeError(Bivector<double, World>(i_),
                    Bivector<double, World>(i_ + j_ + k_)),
      Eq(Sqrt(2)));
}

TEST_F(NumericsTest, TrivectorRelativeError) {
  EXPECT_THAT(
      RelativeError(Trivector<double, World>(42.0),
                    Trivector<double, World>(42.0)),
      Eq(0));
  EXPECT_THAT(
      RelativeError(Trivector<double, World>(1.0),
                    Trivector<double, World>(-1.0)),
      Eq(2));
  EXPECT_THAT(
      RelativeError(Trivector<double, World>(2.0),
                    Trivector<double, World>(1.0)),
      Eq(0.5));
  EXPECT_THAT(
      RelativeError(Trivector<double, World>(1.0),
                    Trivector<double, World>(2.0)),
      Eq(1));
  EXPECT_THAT(
      RelativeError(Trivector<double, World>(42.0),
                    Trivector<double, World>(6.0 * 9.0)),
      AllOf(Gt(0.28), Lt(0.29)));
}

}  // namespace testing_utilities
}  // namespace principia
