#include "geometry/r3x3_matrix.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Eq;

namespace principia {
namespace geometry {

class R3x3MatrixTest : public testing::Test {
 protected:
  using R3 = R3Element<double>;

  void SetUp() override {
    m1_ = R3x3Matrix({-9, 6, 6}, {7, -5, -4}, {-1 ,2, 1});
    m2_ = R3x3Matrix({1, 2, 2}, {-1, -1, 2}, {3, 4, 1});
  }

  R3x3Matrix m1_;
  R3x3Matrix m2_;
  R3x3Matrix m3_;
};

TEST_F(R3x3MatrixTest, Transpose) {
  EXPECT_THAT(m1_.Transpose(),
              Eq(R3x3Matrix({-9, 7, -1}, {6, -5, 2}, {6, -4, 1})));
}

TEST_F(R3x3MatrixTest, UnaryOperators) {
  EXPECT_THAT(+m1_,
              Eq(R3x3Matrix({-9, 6, 6}, {7, -5, -4}, {-1 ,2, 1})));
  EXPECT_THAT(-m2_,
              Eq(R3x3Matrix({-1, -2, -2}, {1, 1, -2}, {-3, -4, -1})));
}

TEST_F(R3x3MatrixTest, BinaryOperators) {
  EXPECT_THAT(m1_ + m2_,
              Eq(R3x3Matrix({-8, 8, 8}, {6, -6, -2}, {2 ,6, 2})));
  EXPECT_THAT(m1_ - m2_,
              Eq(R3x3Matrix({-10, 4, 4}, {8, -4, -6}, {-4, -2, 0})));
  EXPECT_THAT(m1_ * m2_,
              Eq(R3x3Matrix({3, 0, 0}, {0, 3, 0}, {0, 0, 3})));
}

TEST_F(R3x3MatrixTest, ScalarMultiplicationDivision) {
  EXPECT_THAT(3 * m1_,
              Eq(R3x3Matrix({-27, 18, 18}, {21, -15, -12}, {-3 ,6, 3})));
  EXPECT_THAT(m2_ * 5,
              Eq(R3x3Matrix({5, 10, 10}, {-5, -5, 10}, {15, 20, 5})));
  EXPECT_THAT(m1_ / 4,
              Eq(R3x3Matrix({-2.25, 1.5, 1.5},
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
              Eq(R3x3Matrix({-8, 8, 8}, {6, -6, -2}, {2 ,6, 2})));
  EXPECT_THAT(b,
              Eq(R3x3Matrix({-10, 4, 4}, {8, -4, -6}, {-4, -2, 0})));
  EXPECT_THAT(c,
              Eq(R3x3Matrix({3, 0, 0}, {0, 3, 0}, {0, 0, 3})));
  EXPECT_THAT(d,
              Eq(R3x3Matrix({-27, 18, 18}, {21, -15, -12}, {-3 ,6, 3})));
  EXPECT_THAT(e,
              Eq(R3x3Matrix({-2.25, 1.5, 1.5},
                            {1.75, -1.25, -1},
                            {-0.25, 0.5, 0.25})));
}

}  // namespace geometry
}  // namespace principia
