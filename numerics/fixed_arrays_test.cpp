#include "numerics/fixed_arrays.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace numerics {

class FixedArraysTest : public ::testing::Test {
 protected:
  FixedArraysTest()
    : m34_({-8, -6, -4, -7,
            -4, -10, 9, -5,
             6, -3, -2, -9}),
      v3_({10, 31, -47}),
      v4_({-3, -3, 1, 4}),
      l4_({
           1,
           2, 3,
           5, 8, 13}) {}

  FixedMatrix<double, 3, 4> m34_;
  FixedVector<double, 3> v3_;
  FixedVector<double, 4> v4_;
  FixedStrictlyLowerTriangularMatrix<double, 4> l4_;
};

TEST_F(FixedArraysTest, Assignment) {
  FixedVector<double, 2> u2({1, 2});
  FixedVector<double, 2> v2 = {1, 2};
  FixedVector<double, 2> w2;
  w2 = {1, 2};
  EXPECT_EQ(u2, v2);
  EXPECT_EQ(u2, w2);

  FixedMatrix<double, 2, 3> l23({1, 2, 3,
                                4, 5, 6});
  FixedMatrix<double, 2, 3> m23 = {1, 2, 3,
                                   4, 5, 6};
  FixedMatrix<double, 2, 3> n23 = {0, 0, 0,
                                   0, 0, 0};
  n23 = {1, 2, 3,
         4, 5, 6};
  EXPECT_EQ(l23, m23);
  EXPECT_EQ(l23, n23);
}

TEST_F(FixedArraysTest, Multiplication) {
  EXPECT_EQ(v3_, m34_ * v4_);
}

TEST_F(FixedArraysTest, VectorIndexing) {
  EXPECT_EQ(31, v3_[1]);
  v3_[2] = -666;
  EXPECT_EQ(-666, v3_[2]);
}

TEST_F(FixedArraysTest, StrictlyLowerTriangularMatrixIndexing) {
  EXPECT_EQ(6, (FixedStrictlyLowerTriangularMatrix<double, 4>::kDimension));
  EXPECT_EQ(1, l4_[1][0]);
  EXPECT_EQ(2, l4_[2][0]);
  EXPECT_EQ(3, l4_[2][1]);
  EXPECT_EQ(5, l4_[3][0]);
  EXPECT_EQ(8, l4_[3][1]);
  EXPECT_EQ(13, l4_[3][2]);
  l4_[3][1] = -666;
  EXPECT_EQ(-666, l4_[3][1]);
}

}  // namespace numerics
}  // namespace principia
