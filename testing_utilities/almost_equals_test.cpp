#include "testing_utilities/almost_equals.hpp"

#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/bipm.hpp"
#include "quantities/cgs.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/uk.hpp"

namespace principia {
namespace testing_utilities {

using bipm::Knot;
using geometry::Bivector;
using geometry::R3Element;
using geometry::Vector;
using geometry::Trivector;
using quantities::Length;
using quantities::MagneticFlux;
using quantities::Speed;
using testing::Ne;
using testing::Eq;
using testing::Not;
using uk::Foot;
using cgs::Maxwell;

namespace {
struct World;
}  // namespace

class AlmostEqualsTest : public testing::Test {};

TEST_F(AlmostEqualsTest, Dimensionless) {
  double const y = e;
  EXPECT_THAT(y, AlmostEquals(e));
  EXPECT_THAT(y, AlmostEquals(e, 0));
  EXPECT_THAT(2 * y, Not(AlmostEquals(y)));
  double const δy = e / 100.0;
  double e_accumulated = 0.0;
  for (int i = 1; i <= 100.0; ++i) {
    e_accumulated += δy;
  }
  EXPECT_THAT(e_accumulated, Ne(e));
  EXPECT_THAT(e_accumulated, Not(AlmostEquals(e, 0)));
  EXPECT_THAT(e_accumulated, AlmostEquals(e));
}

TEST_F(AlmostEqualsTest, Quantity) {
  Speed v1 = 1 * Knot;
  Speed const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1)));
  Speed const δv = v1 / 100;
  Speed v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1));
}

TEST_F(AlmostEqualsTest, R3Element) {
  R3Element<Speed> const v1 = {1 * Knot, 2 * Knot, 3 * Knot};
  R3Element<Speed> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1)));
  R3Element<Speed> const δv = v1 / 100;
  R3Element<Speed> v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 8));
}

TEST_F(AlmostEqualsTest, Vector) {
  Vector<Length, World> const v1({1 * Foot, 2 * Foot, 3 * Foot});
  Vector<Length, World> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1)));
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
  EXPECT_THAT(v2, AlmostEquals(v1));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1)));
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
  EXPECT_THAT(v2, AlmostEquals(v1));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1)));
  Trivector<MagneticFlux, World> const δv = v1 / 100;
  Trivector<MagneticFlux, World> v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 9));
}

}  // namespace testing_utilities
}  // namespace principia
