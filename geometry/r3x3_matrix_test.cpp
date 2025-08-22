#include "geometry/r3x3_matrix.hpp"

#include <utility>

#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using ::testing::Eq;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

class R3x3MatrixTest : public ::testing::Test {
 protected:
  using R3 = R3Element<double>;

  R3x3MatrixTest()
      : m1_(R3x3Matrix<double>({-9, 6, 6}, {7, -5, -4}, {-1, 2, 1})),
        m2_(R3x3Matrix<double>({1, 2, 2}, {-1, -1, 2}, {3, 4, 1})) {}

  R3x3Matrix<double> m1_;
  R3x3Matrix<double> m2_;
};

using R3x3MatrixDeathTest = R3x3MatrixTest;

TEST_F(R3x3MatrixTest, Trace) {
  EXPECT_THAT(m1_.Trace(), Eq(-13));
}

TEST_F(R3x3MatrixTest, Determinant) {
  EXPECT_THAT(m1_.Determinant(), Eq(9));
  EXPECT_THAT(m2_.Determinant(), Eq(3));
}

TEST_F(R3x3MatrixTest, Transpose) {
  EXPECT_THAT(m1_.Transpose(),
              Eq(R3x3Matrix<double>({-9, 7, -1}, {6, -5, 2}, {6, -4, 1})));
}

TEST_F(R3x3MatrixTest, FrobeniusNorm) {
  EXPECT_THAT(m1_.FrobeniusNorm(), Eq(Sqrt(249)));
}

TEST_F(R3x3MatrixTest, QRDecomposition) {
  R3x3Matrix<Length> hilbert(
      {1 * Metre, 1.0 / 2.0 * Metre, 1.0 / 3.0 * Metre},
      {1.0 / 2.0 * Metre, 1.0 / 3.0 * Metre, 1.0 / 4.0 * Metre},
      {1.0 / 3.0 * Metre, 1.0 / 4.0 * Metre, 1.0 / 5.0 * Metre});
  R3x3Matrix<double> q;
  R3x3Matrix<Length> r;
  hilbert.QRDecomposition(q, r);
  EXPECT_THAT(
      q,
      AlmostEquals(R3x3Matrix<double>(
          {6.0 / 7.0, -30.0 / Sqrt(3577.0),  1.0 / Sqrt(73.0)},
          {3.0 / 7.0,  34.0 / Sqrt(3577.0), -6.0 / Sqrt(73.0)},
          {2.0 / 7.0,  39.0 / Sqrt(3577.0), 6.0 / Sqrt(73.0)}), 569));
  EXPECT_THAT(
      r,
      AlmostEquals(R3x3Matrix<Length>(
          {7.0 / 6.0 * Metre, 9.0 / 14.0 * Metre, 9.0 / 20.0 * Metre},
          {0 * Metre, Sqrt(73.0) / 84.0 * Metre, 9.0 / Sqrt(7300.0) * Metre},
          {0 * Metre, 0 * Metre, 1 / Sqrt(65700.0) * Metre}), 23));
  EXPECT_THAT(q * r, AlmostEquals(hilbert, 0));
}

TEST_F(R3x3MatrixTest, Solve) {
  {
    R3x3Matrix<double> a({2, -3, -4}, {0, 0, -1}, {1, -2, 1});
    R3Element<Length> b(2 * Metre, 5 * Metre, 3 * Metre);
    EXPECT_THAT(a.Solve(b),
                Eq(R3Element<Length>(-60 * Metre, -34 * Metre, -5 * Metre)));
  }
  {
    R3x3Matrix<double> hilbert({1, 1.0 / 2.0, 1.0 / 3.0},
                               {1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0},
                               {1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0});
    R3Element<double> b(6, -12, 5);
    EXPECT_THAT(hilbert.Solve(b),
                AlmostEquals(R3Element<double>(636, -3420, 3240), 48));
  }
  {
    R3x3Matrix<double> vandermonde({1, 2, 4}, {1, -3, 9}, {1, 5, 25});
    R3Element<double> b(7, -9, 11);
    EXPECT_THAT(
        vandermonde.Solve(b),
        AlmostEquals(R3Element<double>(2, 89.0 / 30.0, -7.0 / 30.0), 4));
  }
}

#ifdef _DEBUG
TEST_F(R3x3MatrixDeathTest, IndexingError) {
  EXPECT_DEATH({
    m1_(-1, 2);
  }, "indices = \\{-1, 2\\}");
  EXPECT_DEATH({
    m1_(2, -1);
  }, "index = -1");
  EXPECT_DEATH({
    m1_(1, 3);
  }, "index = 3");
  EXPECT_DEATH({
    m1_(3, 1);
  }, "indices = \\{3, 1\\}");
}
#endif

TEST_F(R3x3MatrixTest, IndexingSuccess) {
  double const a = m1_(1, 2);
  double const b = m1_(0, 0);
  EXPECT_THAT(a, Eq(-4));
  EXPECT_THAT(b, Eq(-9));
}

TEST_F(R3x3MatrixTest, UnaryOperators) {
  EXPECT_THAT(+m1_,
              Eq(R3x3Matrix<double>({-9, 6, 6}, {7, -5, -4}, {-1, 2, 1})));
  EXPECT_THAT(-m2_,
              Eq(R3x3Matrix<double>({-1, -2, -2}, {1, 1, -2}, {-3, -4, -1})));
}

TEST_F(R3x3MatrixTest, BinaryOperators) {
  EXPECT_THAT(m1_ + m2_,
              Eq(R3x3Matrix<double>({-8, 8, 8}, {6, -6, -2}, {2, 6, 2})));
  EXPECT_THAT(m1_ - m2_,
              Eq(R3x3Matrix<double>({-10, 4, 4}, {8, -4, -6}, {-4, -2, 0})));
  EXPECT_THAT(m1_ * m2_,
              Eq(R3x3Matrix<double>({3, 0, 0}, {0, 3, 0}, {0, 0, 3})));
}

TEST_F(R3x3MatrixTest, ScalarMultiplicationDivision) {
  EXPECT_THAT((3 * Metre) * m1_,
              Eq(R3x3Matrix<Length>({-27 * Metre, 18 * Metre, 18 * Metre},
                                    {21 * Metre, -15 * Metre, -12 * Metre},
                                    {-3 * Metre, 6 * Metre, 3 * Metre})));
  EXPECT_THAT(m2_ * 5,
              Eq(R3x3Matrix<double>({5, 10, 10}, {-5, -5, 10}, {15, 20, 5})));
  EXPECT_THAT(m1_ / 4,
              Eq(R3x3Matrix<double>({-2.25, 1.5, 1.5},
                                    {1.75, -1.25, -1},
                                    {-0.25, 0.5, 0.25})));
}

TEST_F(R3x3MatrixTest, Assignment) {
  R3x3Matrix<double> a = m1_;
  R3x3Matrix<double> b = m1_;
  R3x3Matrix<double> c = m1_;
  R3x3Matrix<double> d = m1_;
  R3x3Matrix<double> e = m1_;
  a += m2_;
  b -= m2_;
  c *= m2_;
  d *= 3;
  e /= 4;
  EXPECT_THAT(a,
              Eq(R3x3Matrix<double>({-8, 8, 8}, {6, -6, -2}, {2, 6, 2})));
  EXPECT_THAT(b,
              Eq(R3x3Matrix<double>({-10, 4, 4}, {8, -4, -6}, {-4, -2, 0})));
  EXPECT_THAT(c,
              Eq(R3x3Matrix<double>({3, 0, 0}, {0, 3, 0}, {0, 0, 3})));
  EXPECT_THAT(d,
              Eq(R3x3Matrix<double>({-27, 18, 18},
                                    {21, -15, -12},
                                    {-3, 6, 3})));
  EXPECT_THAT(e,
              Eq(R3x3Matrix<double>({-2.25, 1.5, 1.5},
                                    {1.75, -1.25, -1},
                                    {-0.25, 0.5, 0.25})));
}

TEST_F(R3x3MatrixTest, KroneckerProduct) {
  R3Element<double> v1{1, 2, 4};
  R3Element<double> v2{2, 3, 5};
  EXPECT_THAT(KroneckerProduct(v1, v2),
              Eq(R3x3Matrix<double>({2,  3,  5},
                                    {4,  6, 10},
                                    {8, 12, 20})));
}

TEST_F(R3x3MatrixTest, Serialization) {
  serialization::R3x3Matrix message;
  m1_.WriteToMessage(&message);
  EXPECT_TRUE(message.row_x().x().has_double_());
  EXPECT_EQ(-9.0, message.row_x().x().double_());
  EXPECT_TRUE(message.row_x().y().has_double_());
  EXPECT_EQ(6.0, message.row_x().y().double_());
  EXPECT_TRUE(message.row_x().z().has_double_());
  EXPECT_EQ(6.0, message.row_x().z().double_());
  EXPECT_TRUE(message.row_y().x().has_double_());
  EXPECT_EQ(7.0, message.row_y().x().double_());
  EXPECT_TRUE(message.row_y().y().has_double_());
  EXPECT_EQ(-5.0, message.row_y().y().double_());
  EXPECT_TRUE(message.row_y().z().has_double_());
  EXPECT_EQ(-4.0, message.row_y().z().double_());
  EXPECT_TRUE(message.row_z().x().has_double_());
  EXPECT_EQ(-1.0, message.row_z().x().double_());
  EXPECT_TRUE(message.row_z().y().has_double_());
  EXPECT_EQ(2.0, message.row_z().y().double_());
  EXPECT_TRUE(message.row_z().z().has_double_());
  EXPECT_EQ(1.0, message.row_z().z().double_());
  auto const m = R3x3Matrix<double>::ReadFromMessage(message);
  EXPECT_EQ(m1_, m);
}

}  // namespace geometry
}  // namespace principia
