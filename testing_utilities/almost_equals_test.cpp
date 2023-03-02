#include "testing_utilities/almost_equals.hpp"

#include <sstream>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/bipm.hpp"
#include "quantities/cgs.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"

namespace principia {
namespace testing_utilities {

using numerics::FixedLowerTriangularMatrix;
using numerics::FixedUpperTriangularMatrix;
using numerics::FixedVector;
using numerics::UnboundedLowerTriangularMatrix;
using numerics::UnboundedUpperTriangularMatrix;
using numerics::UnboundedVector;
using testing::Ne;
using testing::Eq;
using testing::Not;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_rotation;
using namespace principia::quantities::_bipm;
using namespace principia::quantities::_cgs;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_uk;
namespace si = quantities::si;

namespace {
using World = Frame<struct WorldTag>;
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

TEST_F(AlmostEqualsTest, R3x3Matrix) {
  R3x3Matrix<Speed> const m1({1 * Knot, 2 * Knot, 3 * Knot},
                             {4 * Knot, -5 * Knot, 7 * Knot},
                             {10 * Knot, 2 * Knot, -30 * Knot});
  R3x3Matrix<Speed> const m2 = m1;
  EXPECT_THAT(m2, AlmostEquals(m1, 0));
  EXPECT_THAT(2 * m2, Not(AlmostEquals(m1, 4)));
  R3x3Matrix<Speed> const δm = m1 / 100;
  R3x3Matrix<Speed> m_accumulated;
  for (int i = 1; i <= 100; ++i) {
    m_accumulated += δm;
  }
  EXPECT_THAT(m_accumulated, Ne(m1));
  EXPECT_THAT(m_accumulated, AlmostEquals(m1, 16));
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

TEST_F(AlmostEqualsTest, Rotation) {
  Rotation<World, World> const r1(Quaternion{1, {2, 3, 4}});
  Rotation<World, World> const r2(Quaternion{-1, {-2, -3, -4}});
  EXPECT_THAT(r2, AlmostEquals(r1, 0));
  EXPECT_THAT(r1 * r2, Not(AlmostEquals(r1, 4)));
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

TEST_F(AlmostEqualsTest, FixedVector) {
  FixedVector<double, 3> const v1({1, 2, 3});
  FixedVector<double, 3> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1, 0));
  EXPECT_THAT(v2, Not(AlmostEquals(v1, 4)));
  double const δv = v1[1] / 100;
  FixedVector<double, 3> v_accumulated({1, 0, 3});
  for (int i = 1; i <= 100; ++i) {
    v_accumulated[1] += δv;
  }
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 3));
}

TEST_F(AlmostEqualsTest, FixedLowerTriangularMatrix) {
  FixedLowerTriangularMatrix<double, 3> const m1({1,
                                               2, 3,
                                               4, 5, 6});
  FixedLowerTriangularMatrix<double, 3> const m2 = m1;
  EXPECT_THAT(m2, AlmostEquals(m1, 0));
  EXPECT_THAT(m2, Not(AlmostEquals(m1, 4)));
  double const δv = m1(1, 0) / 100;
  FixedLowerTriangularMatrix<double, 3> m_accumulated({1,
                                                       0, 3,
                                                       4, 5, 6});
  for (int i = 1; i <= 100; ++i) {
    m_accumulated(1, 0) += δv;
  }
  EXPECT_THAT(m_accumulated, AlmostEquals(m1, 3));
}

TEST_F(AlmostEqualsTest, FixedUpperTriangularMatrix) {
  FixedUpperTriangularMatrix<double, 3> const m1({1, 2, 3,
                                                     4, 5,
                                                        6});
  FixedUpperTriangularMatrix<double, 3> const m2 = m1;
  EXPECT_THAT(m2, AlmostEquals(m1, 0));
  EXPECT_THAT(m2, Not(AlmostEquals(m1, 4)));
  double const δv = m1(0, 1) / 100;
  FixedUpperTriangularMatrix<double, 3> m_accumulated({1, 0, 3,
                                                          4, 5,
                                                             6});
  for (int i = 1; i <= 100; ++i) {
    m_accumulated(0, 1) += δv;
  }
  EXPECT_THAT(m_accumulated, AlmostEquals(m1, 3));
}

TEST_F(AlmostEqualsTest, UnboundedVector) {
  UnboundedVector<double> const v1({1, 2, 3});
  UnboundedVector<double> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1, 0));
  EXPECT_THAT(v2, Not(AlmostEquals(v1, 4)));
  double const δv = v1[1] / 100;
  UnboundedVector<double> v_accumulated({1, 0, 3});
  for (int i = 1; i <= 100; ++i) {
    v_accumulated[1] += δv;
  }
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 3));
}

TEST_F(AlmostEqualsTest, UnboundedLowerTriangularMatrix) {
  UnboundedLowerTriangularMatrix<double> const m1({1,
                                                   2, 3,
                                                   4, 5, 6});
  UnboundedLowerTriangularMatrix<double> const m2 = m1;
  EXPECT_THAT(m2, AlmostEquals(m1, 0));
  EXPECT_THAT(m2, Not(AlmostEquals(m1, 4)));
  double const δv = m1(1, 0) / 100;
  UnboundedLowerTriangularMatrix<double> m_accumulated({1,
                                                        0, 3,
                                                        4, 5, 6});
  for (int i = 1; i <= 100; ++i) {
    m_accumulated(1, 0) += δv;
  }
  EXPECT_THAT(m_accumulated, AlmostEquals(m1, 3));
}

TEST_F(AlmostEqualsTest, UnboundedUpperTriangularMatrix) {
  UnboundedUpperTriangularMatrix<double> const m1({1, 2, 3,
                                                      4, 5,
                                                         6});
  UnboundedUpperTriangularMatrix<double> const m2 = m1;
  EXPECT_THAT(m2, AlmostEquals(m1, 0));
  EXPECT_THAT(m2, Not(AlmostEquals(m1, 4)));
  double const δv = m1(0, 1) / 100;
  UnboundedUpperTriangularMatrix<double> m_accumulated({1, 0, 3,
                                                           4, 5,
                                                              6});
  for (int i = 1; i <= 100; ++i) {
    m_accumulated(0, 1) += δv;
  }
  EXPECT_THAT(m_accumulated, AlmostEquals(m1, 3));
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
