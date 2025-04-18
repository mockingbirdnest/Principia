#include "numerics/fixed_arrays.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/transposed_view.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_elementary_functions;
using ::testing::Pointer;

class FixedArraysTest : public ::testing::Test {
 protected:
  FixedArraysTest()
    : u3_({6, -1, 12}),
      v3_({10, 31, -47}),
      v4_({-3, -3, 1, 4}),
      m34_({-8,  -6, -4, -7,
            -4, -10,  9, -5,
             6,  -3, -2, -9}),
      m23_({1, -2,  0,
            2,  3,  7}),
      n23_({ 5, -1,  3,
            12, 13, -4}),
      sl4_({
            1,
            2, 3,
            5, 8, 13}),
      l4_({ 1,
            2,  3,
            5,  8,  13,
            21, 34, 55, 89}),
      su4_({    1, 2,  3,
                   5,  8,
                      13
                        }),
      u4_({1,  2,  3,  5,
               8, 13, 21,
                  34, 55,
                      89}) {}

  FixedVector<double, 3> u3_;
  FixedVector<double, 3> v3_;
  FixedVector<double, 4> v4_;
  FixedMatrix<double, 3, 4> m34_;
  FixedMatrix<double, 2, 3> m23_;
  FixedMatrix<double, 2, 3> n23_;
  FixedStrictlyLowerTriangularMatrix<double, 4> sl4_;
  FixedLowerTriangularMatrix<double, 4> l4_;
  FixedStrictlyUpperTriangularMatrix<double, 4> su4_;
  FixedUpperTriangularMatrix<double, 4> u4_;
};

TEST_F(FixedArraysTest, Assignment) {
  {
    FixedVector<double, 2> u2({1, 2});
    FixedVector<double, 2> v2 = {{1, 2}};
    FixedVector<double, 2> w2;
    w2 = {1, 2};
    EXPECT_EQ(u2, v2);
    EXPECT_EQ(u2, w2);
  }

  {
    FixedMatrix<double, 2, 3> l23({1, 2, 3,
                                   4, 5, 6});
    FixedMatrix<double, 2, 3> m23 = {{1, 2, 3,
                                      4, 5, 6}};
    FixedMatrix<double, 2, 3> n23 = {{0, 0, 0,
                                      0, 0, 0}};
    n23 = {1, 2, 3,
           4, 5, 6};
    EXPECT_EQ(l23, m23);
    EXPECT_EQ(l23, n23);
  }
  {
    FixedStrictlyLowerTriangularMatrix<double, 3> l3({
                                                      1,
                                                      2, 3 });
    FixedStrictlyLowerTriangularMatrix<double, 3> m3 = {{
                                                         1,
                                                         2, 3}};
    FixedStrictlyLowerTriangularMatrix<double, 3> n3 = {{
                                                         0,
                                                         0, 0}};
    FixedStrictlyLowerTriangularMatrix<double, 3> o3;
    EXPECT_EQ(o3, n3);
    n3 = {
          1,
          2, 3};
    EXPECT_EQ(l3, m3);
    EXPECT_EQ(l3, n3);
  }
  {
    FixedLowerTriangularMatrix<double, 3> l3({1,
                                              2, 3,
                                              4, 5, 6});
    FixedLowerTriangularMatrix<double, 3> m3 = {{1,
                                                 2, 3,
                                                 4, 5, 6}};
    FixedLowerTriangularMatrix<double, 3> n3 = {{0,
                                                 0, 0,
                                                 0, 0, 0}};
    FixedLowerTriangularMatrix<double, 3> o3;
    EXPECT_EQ(o3, n3);
    n3 = {1,
          2, 3,
          4, 5, 6};
    EXPECT_EQ(l3, m3);
    EXPECT_EQ(l3, n3);
  }
  {
    FixedUpperTriangularMatrix<double, 3> l3({1, 2, 3,
                                                 4, 5,
                                                    6});
    FixedUpperTriangularMatrix<double, 3> m3 = {{1, 2, 3,
                                                    4, 5,
                                                       6}};
    FixedUpperTriangularMatrix<double, 3> n3 = {{0, 0, 0,
                                                    0, 0,
                                                       0}};
    FixedUpperTriangularMatrix<double, 3> o3;
    EXPECT_EQ(o3, n3);
    n3 = {1, 2, 3,
             4, 5,
                6};
    EXPECT_EQ(l3, m3);
    EXPECT_EQ(l3, n3);
  }
}

TEST_F(FixedArraysTest, Norm) {
  EXPECT_EQ(35, TransposedView{v4_} * v4_);  // NOLINT
  EXPECT_EQ(Sqrt(35.0), v4_.Norm());
  EXPECT_EQ(35, v4_.Norm²());
  EXPECT_EQ(Sqrt(517.0), m34_.FrobeniusNorm());
}

TEST_F(FixedArraysTest, AdditiveGroups) {
  EXPECT_EQ((FixedVector<double, 3>({-10, -31, 47})), -v3_);
  EXPECT_EQ((FixedMatrix<double, 2, 3>({-1,  2,  0,
                                        -2, -3, -7})), -m23_);

  EXPECT_EQ((FixedVector<double, 3>({16, 30, -35})), u3_ + v3_);
  EXPECT_EQ((FixedMatrix<double, 2, 3>({ 6, -3, 3,
                                        14, 16, 3})), m23_ + n23_);

  EXPECT_EQ((FixedVector<double, 3>({-4, -32, 59})), u3_ - v3_);
  EXPECT_EQ((FixedMatrix<double, 2, 3>({ -4,  -1, -3,
                                        -10, -10, 11})), m23_ - n23_);
}

TEST_F(FixedArraysTest, VectorSpaces) {
  EXPECT_EQ((FixedVector<double, 3>({12, -2, 24})), 2 * u3_);
  EXPECT_EQ((FixedVector<double, 3>({-30, -93, 141})), v3_ * -3);

  EXPECT_EQ((FixedMatrix<double, 2, 3>({2, -4,  0,
                                        4,  6, 14})), 2 * m23_);
  EXPECT_EQ((FixedMatrix<double, 2, 3>({-15,   3, -9,
                                        -36, -39, 12})), n23_ * -3);

  EXPECT_EQ((FixedVector<double, 3>({2.5, 7.75, -11.75})), v3_ / 4);
  EXPECT_EQ((FixedMatrix<double, 2, 3>({-2.5,  0.5, -1.5,
                                          -6, -6.5,  2})), n23_ / -2);
}

TEST_F(FixedArraysTest, Algebra) {
  EXPECT_EQ(-535, TransposedView{u3_} * v3_);  // NOLINT
  EXPECT_EQ((FixedMatrix<double, 3, 4>({-30, -30,  10,   40,
                                        -93, -93,  31,  124,
                                        141, 141, -47, -188})),
             v3_ * TransposedView{v4_});
  EXPECT_EQ((FixedMatrix<double, 2, 4>({ 0,  14, -22,   3,
                                        14, -63,   5, -92})),
            m23_ * m34_);
  auto m = FixedMatrix<double, 2, 2>({1, 2,
                                      3, 4});
  EXPECT_THAT(&(m *= FixedMatrix<double, 2, 2>({5, 6,
                                                7, 8})),
              Pointer(&m));
  EXPECT_EQ((FixedMatrix<double, 2, 2>({19, 22,
                                        43, 50})),
            m);
  EXPECT_EQ(v3_, m34_ * v4_);
  EXPECT_EQ((FixedVector<double, 4>({-486, -229, 333, 198})),
            TransposedView{m34_} * v3_);  // NOLINT
  FixedMatrix<double, 4, 3> m43(TransposedView{m34_});
  EXPECT_EQ((FixedVector<double, 4>({-486, -229, 333, 198})),
            m43 * v3_);
}

TEST_F(FixedArraysTest, VectorIndexing) {
  EXPECT_EQ(31, v3_[1]);
  v3_[2] = -666;
  EXPECT_EQ(-666, v3_[2]);
}

TEST_F(FixedArraysTest, MatrixIndexing) {
  EXPECT_EQ(9, m34_(1, 2));
  m34_(2, 1) = -666;
  EXPECT_EQ(-666, m34_(2, 1));
}

TEST_F(FixedArraysTest, StrictlyLowerTriangularMatrixIndexing) {
  EXPECT_EQ(1, sl4_(1, 0));
  EXPECT_EQ(2, sl4_(2, 0));
  EXPECT_EQ(3, sl4_(2, 1));
  EXPECT_EQ(5, sl4_(3, 0));
  EXPECT_EQ(8, sl4_(3, 1));
  EXPECT_EQ(13, sl4_(3, 2));
  sl4_(3, 1) = -666;
  EXPECT_EQ(-666, sl4_(3, 1));

  FixedStrictlyLowerTriangularMatrix<double, 4> const sl4 = sl4_;
  EXPECT_EQ(1, sl4(1, 0));
}

TEST_F(FixedArraysTest, LowerTriangularMatrixIndexing) {
  EXPECT_EQ(1, l4_(0, 0));
  EXPECT_EQ(2, l4_(1, 0));
  EXPECT_EQ(3, l4_(1, 1));
  EXPECT_EQ(5, l4_(2, 0));
  EXPECT_EQ(8, l4_(2, 1));
  EXPECT_EQ(13, l4_(2, 2));
  EXPECT_EQ(21, l4_(3, 0));
  EXPECT_EQ(34, l4_(3, 1));
  EXPECT_EQ(55, l4_(3, 2));
  EXPECT_EQ(89, l4_(3, 3));
  l4_(3, 1) = -666;
  EXPECT_EQ(-666, l4_(3, 1));

  FixedLowerTriangularMatrix<double, 4> const l4 = l4_;
  EXPECT_EQ(1, l4(0, 0));
}

TEST_F(FixedArraysTest, StrictlyUpperTriangularMatrixIndexing) {
  EXPECT_EQ(1, su4_(0, 1));
  EXPECT_EQ(2, su4_(0, 2));
  EXPECT_EQ(3, su4_(0, 3));
  EXPECT_EQ(5, su4_(1, 2));
  EXPECT_EQ(8, su4_(1, 3));
  EXPECT_EQ(13, su4_(2, 3));
  su4_(1, 3) = -666;
  EXPECT_EQ(-666, su4_(1, 3));

  FixedStrictlyUpperTriangularMatrix<double, 4> const su4 = su4_;
  EXPECT_EQ(1, su4(0, 1));
}

TEST_F(FixedArraysTest, UpperTriangularMatrixIndexing) {
  EXPECT_EQ(1, u4_(0, 0));
  EXPECT_EQ(2, u4_(0, 1));
  EXPECT_EQ(3, u4_(0, 2));
  EXPECT_EQ(5, u4_(0, 3));
  EXPECT_EQ(8, u4_(1, 1));
  EXPECT_EQ(13, u4_(1, 2));
  EXPECT_EQ(21, u4_(1, 3));
  EXPECT_EQ(34, u4_(2, 2));
  EXPECT_EQ(55, u4_(2, 3));
  EXPECT_EQ(89, u4_(3, 3));
  u4_(1, 3) = -666;
  EXPECT_EQ(-666, u4_(1, 3));

  FixedUpperTriangularMatrix<double, 4> const u4 = u4_;
  EXPECT_EQ(1, u4(0, 0));
}

TEST_F(FixedArraysTest, Row) {
  FixedMatrix<double, 2, 3> m({1, 2, 3,
                               4, -5, 6});
  auto const* const r0 = m.row<0>();
  auto const* const r1 = m.row<1>();
  FixedVector<double, 3> v = {{1, 2, -3}};

  EXPECT_EQ(-4, r0 * v);
  EXPECT_EQ(-24, r1 * v);
}

TEST_F(FixedArraysTest, Conversions) {
  using M = FixedMatrix<double, 4, 4>;
  EXPECT_EQ(M({0, 0,  0, 0,
               1, 0,  0, 0,
               2, 3,  0, 0,
               5, 8, 13, 0}),
            M(sl4_));
  EXPECT_EQ(M({1,   0,  0,  0,
               2,   3,  0,  0,
               5,   8, 13,  0,
               21, 34, 55, 89}),
            M(l4_));
  EXPECT_EQ(M({0, 1, 2,  3,
               0, 0, 5,  8,
               0, 0, 0, 13,
               0, 0, 0,  0}),
            M(su4_));
  EXPECT_EQ(M({1, 2,  3,  5,
               0, 8, 13, 21,
               0, 0, 34, 55,
               0, 0,  0, 89}),
            M(u4_));
}

TEST_F(FixedArraysTest, Transpose) {
  EXPECT_EQ(
      (FixedMatrix<double, 4, 3>({-8,  -4,  6,
                                  -6, -10, -3,
                                  -4,   9, -2,
                                  -7,  -5, -9})),
      (FixedMatrix<double, 4, 3>(TransposedView{m34_})));
  EXPECT_EQ(
      (FixedUpperTriangularMatrix<double, 4>({1, 2,  5, 21,
                                                 3,  8, 34,
                                                    13, 55,
                                                        89})),
      (FixedUpperTriangularMatrix<double, 4>(TransposedView{l4_})));
  EXPECT_EQ(
      (FixedLowerTriangularMatrix<double, 4>({1,
                                              2,  8,
                                              3, 13, 34,
                                              5, 21, 55, 89})),
      (FixedLowerTriangularMatrix<double, 4>(TransposedView{u4_})));
}

}  // namespace numerics
}  // namespace principia
