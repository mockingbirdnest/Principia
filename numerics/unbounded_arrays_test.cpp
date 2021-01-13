
#include "numerics/unbounded_arrays.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using quantities::Sqrt;
using testing_utilities::AlmostEquals;

class UnboundedArraysTest : public ::testing::Test {
 protected:
  UnboundedArraysTest()
    : v3_({10, 31, -47}),
      v4_({-3, -3, 1, 4}),
      l4_({ 1,
            2,  3,
            5,  8,  13,
            21, 34, 55, 89}),
      u4_({1, 2,  3,  5,
              8, 13, 21,
                 34, 55,
                     89}) {}

  UnboundedVector<double> v3_;
  UnboundedVector<double> v4_;
  UnboundedLowerTriangularMatrix<double> l4_;
  UnboundedUpperTriangularMatrix<double> u4_;
};

TEST_F(UnboundedArraysTest, Assignment) {
  {
    UnboundedVector<double> u2({1, 2});
    UnboundedVector<double> v2 = {{1, 2}};
    UnboundedVector<double> w2(2);
    w2 = {{1, 2}};
    EXPECT_EQ(u2, v2);
    EXPECT_EQ(u2, w2);
  }
  {
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
  {
    UnboundedUpperTriangularMatrix<double> l3({1, 2, 3,
                                                  4, 5,
                                                     6});
    UnboundedUpperTriangularMatrix<double> m3 = {{1, 2, 3,
                                                     4, 5,
                                                        6}};
    UnboundedUpperTriangularMatrix<double> n3 = {{0, 0, 0,
                                                     0, 0,
                                                        0}};
    UnboundedUpperTriangularMatrix<double> o3(3);
    EXPECT_EQ(o3, n3);
    n3 = {{1, 2, 3,
              4, 5,
                 6}};
    EXPECT_EQ(l3, m3);
    EXPECT_EQ(l3, n3);
  }
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

  UnboundedLowerTriangularMatrix<double> const l4 = l4_;
  EXPECT_EQ(1, l4[0][0]);
}

TEST_F(UnboundedArraysTest, UpperTriangularMatrixIndexing) {
  EXPECT_EQ(10, u4_.dimension());
  EXPECT_EQ(1, u4_[0][0]);
  EXPECT_EQ(2, u4_[0][1]);
  EXPECT_EQ(3, u4_[0][2]);
  EXPECT_EQ(5, u4_[0][3]);
  EXPECT_EQ(8, u4_[1][1]);
  EXPECT_EQ(13, u4_[1][2]);
  EXPECT_EQ(21, u4_[1][3]);
  EXPECT_EQ(34, u4_[2][2]);
  EXPECT_EQ(55, u4_[2][3]);
  EXPECT_EQ(89, u4_[3][3]);
  u4_[1][3] = -666;
  EXPECT_EQ(-666, u4_[1][3]);

  UnboundedUpperTriangularMatrix<double> const u4 = u4_;
  EXPECT_EQ(1, u4[0][0]);
}

TEST_F(UnboundedArraysTest, Extend) {
  {
    UnboundedVector<double> u2({1, 2});
    UnboundedVector<double> u4({1, 2, 3, 4});
    u2.Extend({3, 4});
    EXPECT_EQ(u2, u4);
  }
  {
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
  {
    UnboundedUpperTriangularMatrix<double> u3({1, 2, 3,
                                                  4, 5,
                                                     6});
    UnboundedUpperTriangularMatrix<double> u5({1, 2, 3,  7,  8,
                                                  4, 5,  9, 10,
                                                     6, 11, 12,
                                                        13, 14,
                                                            15});
    u3.Extend({ 7,  8,
                9, 10,
               11, 12,
               13, 14,
                   15});
    EXPECT_EQ(u3, u5);
  }
}

TEST_F(UnboundedArraysTest, Erase) {
  {
    UnboundedVector<double> u4({1, 2, 3, 4});
    UnboundedVector<double> u2({1, 2});
    u4.EraseToEnd(2);
    EXPECT_EQ(u2, u4);
  }
  {
    UnboundedLowerTriangularMatrix<double> l4({1,
                                               2, 3,
                                               4, 5, 6,
                                               7, 8, 9, 10});
    UnboundedLowerTriangularMatrix<double> l2({1,
                                               2, 3});
    l4.EraseToEnd(2);
    EXPECT_EQ(l2, l4);
  }
  {
    UnboundedUpperTriangularMatrix<double> u4({1, 2, 3, 4,
                                                  5, 6, 7,
                                                     8, 9,
                                                       10});
    UnboundedUpperTriangularMatrix<double> u2({1, 2,
                                                  5});
    u4.EraseToEnd(2);
    EXPECT_EQ(u2, u4);
  }
}

TEST_F(UnboundedArraysTest, CholeskyFactorization) {
  UnboundedUpperTriangularMatrix<double> const hilbert4({
      1, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0,
         1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0,
                    1.0 / 5.0, 1.0 / 6.0,
                               1.0 / 7.0});
  auto const r4_actual = CholeskyFactorization(hilbert4);
  UnboundedUpperTriangularMatrix<double> const r4_expected({
      1,        1.0 / 2.0,         1.0 / 3.0,          1.0 / 4.0,
         1.0 / Sqrt(12.0),  1.0 / Sqrt(12.0),  Sqrt(27.0) / 20.0,
                           1.0 / Sqrt(180.0),   1.0 / Sqrt(80.0),
                                              1.0 / Sqrt(2800.0)});
  EXPECT_THAT(r4_actual, AlmostEquals(r4_expected, 245));
}

TEST_F(UnboundedArraysTest, ᵗRDRDecomposition) {
  UnboundedUpperTriangularMatrix<double> const hilbert4({
      1, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0,
         1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0,
                    1.0 / 5.0, 1.0 / 6.0,
                               1.0 / 7.0});
  UnboundedUpperTriangularMatrix<double> r4_actual(4);
  UnboundedVector<double> d4_actual(4);
  ᵗRDRDecomposition(hilbert4, r4_actual, d4_actual);
  UnboundedUpperTriangularMatrix<double> const r4_expected({
      1, 1.0 / 2.0, 1.0 / 3.0,  1.0 / 4.0,
                 1,         1, 9.0 / 10.0,
                            1,  3.0 / 2.0,
                                        1});
  UnboundedVector<double> d4_expected(
      {1, 1.0 / 12.0, 1.0 / 180.0, 1.0 / 2800.0});
  EXPECT_THAT(d4_actual, AlmostEquals(d4_expected, 1615));
  EXPECT_THAT(r4_actual, AlmostEquals(r4_expected, 23));
}

}  // namespace numerics
}  // namespace principia
