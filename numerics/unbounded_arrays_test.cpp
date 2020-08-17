
#include "numerics/unbounded_arrays.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace numerics {

class UnboundedArraysTest : public ::testing::Test {
 protected:
  UnboundedArraysTest()
    : v3_({10, 31, -47}),
      v4_({-3, -3, 1, 4}),
      l4_({ 1,
            2,  3,
            5,  8,  13,
            21, 34, 55, 89}) {}

  UnboundedVector<double> v3_;
  UnboundedVector<double> v4_;
  UnboundedLowerTriangularMatrix<double> l4_;
};

TEST_F(UnboundedArraysTest, Assignment) {
  UnboundedVector<double> u2({1, 2});
  UnboundedVector<double> v2 = {{1, 2}};
  UnboundedVector<double> w2(2);
  w2 = {{1, 2}};
  EXPECT_EQ(u2, v2);
  EXPECT_EQ(u2, w2);

  UnboundedLowerTriangularMatrix<double> l3({1,
                                             2, 3,
                                             4, 5, 6});
  UnboundedLowerTriangularMatrix<double> m3 = {{1,
                                                2, 3,
                                                4, 5, 6}};
  UnboundedLowerTriangularMatrix<double> n3 = {{0,
                                                0, 0,
                                                0, 0, 0}};
  UnboundedLowerTriangularMatrix<double> o3(3);
  EXPECT_EQ(o3, n3);
  n3 = {{1,
         2, 3,
         4, 5, 6}};
  EXPECT_EQ(l3, m3);
  EXPECT_EQ(l3, n3);
}

TEST_F(UnboundedArraysTest, VectorIndexing) {
  EXPECT_EQ(31, v3_[1]);
  v3_[2] = -666;
  EXPECT_EQ(-666, v3_[2]);
}

TEST_F(UnboundedArraysTest, LowerTriangularMatrixIndexing) {
  EXPECT_EQ(10, l4_.dimension());
  EXPECT_EQ(1, l4_[0][0]);
  EXPECT_EQ(2, l4_[1][0]);
  EXPECT_EQ(3, l4_[1][1]);
  EXPECT_EQ(5, l4_[2][0]);
  EXPECT_EQ(8, l4_[2][1]);
  EXPECT_EQ(13, l4_[2][2]);
  EXPECT_EQ(21, l4_[3][0]);
  EXPECT_EQ(34, l4_[3][1]);
  EXPECT_EQ(55, l4_[3][2]);
  EXPECT_EQ(89, l4_[3][3]);
  l4_[3][1] = -666;
  EXPECT_EQ(-666, l4_[3][1]);
}

TEST_F(UnboundedArraysTest, Extend) {
  UnboundedVector<double> u2({1, 2});
  UnboundedVector<double> u4({1, 2, 3, 4});
  u2.Extend({3, 4});
  EXPECT_EQ(u2, u4);

  UnboundedLowerTriangularMatrix<double> l3({1,
                                             2, 3,
                                             4, 5, 6});
  UnboundedLowerTriangularMatrix<double> l4({1,
                                             2, 3,
                                             4, 5, 6,
                                             7, 8, 9, 10});
  l3.Extend({7, 8, 9, 10});
  EXPECT_EQ(l3, l4);
}

TEST_F(UnboundedArraysTest, Erase) {
  UnboundedVector<double> u4({1, 2, 3, 4});
  UnboundedVector<double> u2({1, 2});
  u4.EraseToEnd(2);
  EXPECT_EQ(u2, u4);

  UnboundedLowerTriangularMatrix<double> l4({1,
                                             2, 3,
                                             4, 5, 6,
                                             7, 8, 9, 10});
  UnboundedLowerTriangularMatrix<double> l2({1,
                                             2, 3});
  l4.EraseToEnd(2);
  EXPECT_EQ(l2, l4);
}

}  // namespace numerics
}  // namespace principia
