
#include "testing_utilities/numerics.hpp"

#include <cmath>
#include <limits>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace testing_utilities {

using geometry::Bivector;
using geometry::Frame;
using geometry::Point;
using geometry::R3Element;
using geometry::Trivector;
using geometry::Vector;
using numerics::ULPDistance;
using quantities::Sqrt;
using quantities::si::Metre;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Ne;

class NumericsTest : public testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  R3Element<double> const i_ = {1, 0, 0};
  R3Element<double> const j_ = {0, 1, 0};
  R3Element<double> const k_ = {0, 0, 1};
};

double DoubleAbs(const double x) {
  return std::abs(x);
}

TEST_F(NumericsTest, ULPs) {
  EXPECT_THAT(ULPDistance(1, 1), Eq(0));
  EXPECT_THAT(ULPDistance(+0.0, +0.0), Eq(0));
  EXPECT_THAT(ULPDistance(+0.0, -0.0), Eq(0));
  // |std::numeric_limits<double>::min()| is the smallest positive normalized
  // number.  52 bits of mantissa stand between it and 0, in the form of
  // denormals.
  EXPECT_THAT(ULPDistance(+0.0, std::numeric_limits<double>::min()),
              Eq(std::pow(2, std::numeric_limits<double>::digits - 1)));
  EXPECT_THAT(ULPDistance(-0.0, std::numeric_limits<double>::min()),
              Eq(std::pow(2, std::numeric_limits<double>::digits - 1)));
  EXPECT_THAT(ULPDistance(+0.0, -std::numeric_limits<double>::denorm_min()),
              Eq(1));
  EXPECT_THAT(ULPDistance(-0.0, -std::numeric_limits<double>::denorm_min()),
              Eq(1));
  EXPECT_THAT(ULPDistance(-2, 2), Gt(0));
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

TEST_F(NumericsTest, PointAbsoluteError) {
  Point<Vector<double, World>> p1;
  Point<Vector<double, World>> p2;
  EXPECT_THAT(AbsoluteError(p1, p2 + Vector<double, World>(j_)), Eq(1));
}

TEST_F(NumericsTest, DoubleRelativeError) {
  EXPECT_THAT(RelativeError(42.0, 42.0, DoubleAbs), Eq(0));
  EXPECT_THAT(RelativeError(1.0, -1.0, DoubleAbs), Eq(2));
  EXPECT_THAT(RelativeError(2.0, 1.0, DoubleAbs), Eq(0.5));
  EXPECT_THAT(RelativeError(1.0, 2.0, DoubleAbs), Eq(1));
  EXPECT_THAT(RelativeError(42.0, 6.0 * 9.0, DoubleAbs), IsNear(0.28_⑴));
}

TEST_F(NumericsTest, DimensionlessRelativeError) {
  EXPECT_THAT(RelativeError(42.0, 42.0), Eq(0));
  EXPECT_THAT(RelativeError(1.0, -1.0), Eq(2));
  EXPECT_THAT(RelativeError(2.0, 1.0), Eq(0.5));
  EXPECT_THAT(RelativeError(1.0, 2.0), Eq(1));
  EXPECT_THAT(RelativeError(42.0, 6.0 * 9.0), IsNear(0.28_⑴));
}

TEST_F(NumericsTest, DimensionfulRelativeError) {
  EXPECT_THAT(RelativeError(42 * Metre, 42 * Metre), Eq(0));
  EXPECT_THAT(RelativeError(1 * Metre, -1 * Metre), Eq(2));
  EXPECT_THAT(RelativeError(2 * Metre, 1 * Metre), Eq(0.5));
  EXPECT_THAT(RelativeError(1 * Metre, 2 * Metre), Eq(1));
  EXPECT_THAT(RelativeError(42 * Metre, 6 * 9 * Metre), IsNear(0.28_⑴));
}

TEST_F(NumericsTest, R3ElementRelativeError) {
  EXPECT_THAT(RelativeError(i_ + j_ + k_, i_ + j_ + k_), Eq(0));
  EXPECT_THAT(RelativeError(i_ + j_, i_ + j_ + k_), AlmostEquals(Sqrt(0.5), 1));
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
      AlmostEquals(Sqrt(0.5), 1));
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
      AlmostEquals(Sqrt(0.5), 1));
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
      IsNear(0.28_⑴));
}

}  // namespace testing_utilities
}  // namespace principia
