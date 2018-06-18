
#include <utility>

#include "geometry/r3x3_matrix.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace geometry {
namespace internal_r3x3_matrix {

using testing::Eq;

class R3x3MatrixTest : public testing::Test {
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

TEST_F(R3x3MatrixTest, Transpose) {
  EXPECT_THAT(m1_.Transpose(),
              Eq(R3x3Matrix<double>({-9, 7, -1}, {6, -5, 2}, {6, -4, 1})));
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
  EXPECT_THAT(3 * m1_,
              Eq(R3x3Matrix<double>({-27, 18, 18},
                                    {21, -15, -12},
                                    {-3, 6, 3})));
  EXPECT_THAT(m2_ * 5,
              Eq(R3x3Matrix<double>({5, 10, 10}, {-5, -5, 10}, {15, 20, 5})));
  EXPECT_THAT(m1_ / 4,
              Eq(R3x3Matrix<double>({-2.25, 1.5, 1.5},
                                    {1.75, -1.25, -1},
                                    {-0.25, 0.5, 0.25})));
}

TEST_F(R3x3MatrixTest, Assignment) {
  R3x3Matrix a = m1_;
  R3x3Matrix b = m1_;
  R3x3Matrix c = m1_;
  R3x3Matrix d = m1_;
  R3x3Matrix e = m1_;
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

}  // namespace internal_r3x3_matrix
}  // namespace geometry
}  // namespace principia
