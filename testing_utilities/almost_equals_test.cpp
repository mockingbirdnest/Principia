
#include "testing_utilities/almost_equals.hpp"

#include <sstream>

#include "geometry/grassmann.hpp"
#include "geometry/quaternion.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/bipm.hpp"
#include "quantities/cgs.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"

namespace principia {
namespace testing_utilities {

using geometry::Bivector;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::Vector;
using geometry::Trivector;
using quantities::Length;
using quantities::MagneticFlux;
using quantities::Speed;
using quantities::bipm::Knot;
using quantities::cgs::Maxwell;
using quantities::uk::Foot;
using testing::Ne;
using testing::Eq;
using testing::Not;
namespace si = quantities::si;

namespace {
struct World;
}  // namespace

class AlmostEqualsTest : public testing::Test {};

TEST_F(AlmostEqualsTest, Dimensionless) {
  double const y = e;
  EXPECT_THAT(y, AlmostEquals(e, 0));
  EXPECT_THAT(y, Not(AlmostEquals(e, 1)));
  EXPECT_THAT(2 * y, Not(AlmostEquals(y, 4)));
  double const δy = e / 100.0;
  double e_accumulated = 0.0;
  for (int i = 1; i <= 100.0; ++i) {
    e_accumulated += δy;
  }
  EXPECT_THAT(e_accumulated, Ne(e));
  EXPECT_THAT(e_accumulated, Not(AlmostEquals(e, 0)));
  EXPECT_THAT(e_accumulated, AlmostEquals(e, 1));
}

TEST_F(AlmostEqualsTest, Quantity) {
  Speed v1 = 1 * Knot;
  Speed const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1, 0));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1, 4)));
  Speed const δv = v1 / 100.0;
  Speed v_accumulated;
  for (int i = 1; i <= 100.0; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 4));
}

TEST_F(AlmostEqualsTest, R3Element) {
  R3Element<Speed> const v1 = {1 * Knot, 2 * Knot, 3 * Knot};
  R3Element<Speed> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1, 0));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1, 4)));
  R3Element<Speed> const δv = v1 / 100;
  R3Element<Speed> v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 8));
}

TEST_F(AlmostEqualsTest, Quaternion) {
  Quaternion const q1 = {1, {2, 3, 4}};
  Quaternion const q2 = q1;
  EXPECT_THAT(q2, AlmostEquals(q1, 0));
  EXPECT_THAT(2 * q2, Not(AlmostEquals(q1, 4)));
  Quaternion const δq = q1 / 100;
  Quaternion q_accumulated;
  for (int i = 1; i <= 100; ++i) {
    q_accumulated += δq;
  }
  EXPECT_THAT(q_accumulated, Ne(q1));
  EXPECT_THAT(q_accumulated, AlmostEquals(q1, 11));
}

TEST_F(AlmostEqualsTest, Vector) {
  Vector<Length, World> const v1({1 * Foot, 2 * Foot, 3 * Foot});
  Vector<Length, World> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1, 0));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1, 4)));
  Vector<Length, World> const δv = v1 / 100;
  Vector<Length, World> v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 14));
}

TEST_F(AlmostEqualsTest, Bivector) {
  Bivector<double, World> const v1({4, -5, 6});
  Bivector<double, World> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1, 0));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1, 4)));
  Bivector<double, World> const δv = v1 / 100;
  Bivector<double, World> v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 11));
}

TEST_F(AlmostEqualsTest, Trivector) {
  Trivector<MagneticFlux, World> const v1(2 * Maxwell);
  Trivector<MagneticFlux, World> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1, 0));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1, 4)));
  Trivector<MagneticFlux, World> const δv = v1 / 100;
  Trivector<MagneticFlux, World> v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 9));
}

TEST_F(AlmostEqualsTest, Describe) {
  Speed v1 = 1 * si::Unit<Speed>;
  {
    std::ostringstream out;
    AlmostEquals(v1, 2, 6).impl().DescribeTo(&out);
    EXPECT_EQ("is within 2 to 6 ULPs of +1.00000000000000000e+00 m s^-1",
              out.str());
  }
  {
    std::ostringstream out;
    AlmostEquals(v1, 2, 6).impl().DescribeNegationTo(&out);
    EXPECT_EQ("is not within 2 to 6 ULPs of +1.00000000000000000e+00 m s^-1",
              out.str());
  }
}

}  // namespace testing_utilities
}  // namespace principia
