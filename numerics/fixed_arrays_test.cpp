#include "numerics/fixed_arrays.hpp"

#include "gtest/gtest.h"
#include "numerics/transposed_view.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_elementary_functions;

class FixedArraysTest : public ::testing::Test {
 protected:
  FixedArraysTest()
    : m34_({-8,  -6, -4, -7,
            -4, -10,  9, -5,
             6,  -3, -2, -9}),
      m23_({1, -2,  0,
            2,  3,  7}),
      n23_({5,  -1,  3,
            12, 13, -4}),
      u3_({6, -1, 12}),
      v3_({10, 31, -47}),
      v4_({-3, -3, 1, 4}),
      sl4_({
            1,
            2, 3,
            5, 8, 13}),
      l4_({ 1,
            2,  3,
            5,  8,  13,
            21, 34, 55, 89}),
      u4_({1,  2,  3,  5,
               8, 13, 21,
                  34, 55,
                      89}) {}

  FixedMatrix<double, 3, 4> m34_;
  FixedMatrix<double, 2, 3> m23_;
  FixedMatrix<double, 2, 3> n23_;
  FixedVector<double, 3> u3_;
  FixedVector<double, 3> v3_;
  FixedVector<double, 4> v4_;
  FixedStrictlyLowerTriangularMatrix<double, 4> sl4_;
  FixedLowerTriangularMatrix<double, 4> l4_;
  FixedUpperTriangularMatrix<double, 4> u4_;
};

TEST_F(FixedArraysTest, Assignment) {
  FixedVector<double, 2> u2({1, 2});
  FixedVector<double, 2> v2 = {{1, 2}};
  FixedVector<double, 2> w2;
  w2 = {{1, 2}};
  EXPECT_EQ(u2, v2);
  EXPECT_EQ(u2, w2);

  FixedMatrix<double, 2, 3> l23({1, 2, 3,
                                 4, 5, 6});
  FixedMatrix<double, 2, 3> m23 = {{1, 2, 3,
                                    4, 5, 6}};
  FixedMatrix<double, 2, 3> n23 = {{0, 0, 0,
                                    0, 0, 0}};
  n23 = {{1, 2, 3,
          4, 5, 6}};
  EXPECT_EQ(l23, m23);
  EXPECT_EQ(l23, n23);

  FixedStrictlyLowerTriangularMatrix<double, 3> l3({
                                                    1,
                                                    2, 3});
  FixedStrictlyLowerTriangularMatrix<double, 3> m3 = {{
                                                       1,
                                                       2, 3}};
  FixedStrictlyLowerTriangularMatrix<double, 3> n3 = {{
                                                       0,
                                                       0, 0}};
  FixedStrictlyLowerTriangularMatrix<double, 3> o3;
  EXPECT_EQ(o3, n3);
  n3 = {{
         1,
         2, 3}};
  EXPECT_EQ(l3, m3);
  EXPECT_EQ(l3, n3);
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
  EXPECT_EQ(v3_, m34_ * v4_);
  EXPECT_EQ((FixedVector<double, 4>({-486, -229, 333, 198})),
            TransposedView{m34_} * v3_);  // NOLINT
}

TEST_F(FixedArraysTest, VectorIndexing) {
  EXPECT_EQ(31, v3_[1]);
  v3_[2] = -666;
  EXPECT_EQ(-666, v3_[2]);
}

TEST_F(FixedArraysTest, StrictlyLowerTriangularMatrixIndexing) {
  EXPECT_EQ(6, (FixedStrictlyLowerTriangularMatrix<double, 4>::size()));
  EXPECT_EQ(1, sl4_(1, 0));
  EXPECT_EQ(2, sl4_(2, 0));
  EXPECT_EQ(3, sl4_(2, 1));
  EXPECT_EQ(5, sl4_(3, 0));
  EXPECT_EQ(8, sl4_(3, 1));
  EXPECT_EQ(13, sl4_(3, 2));
  sl4_(3, 1) = -666;
  EXPECT_EQ(-666, sl4_(3, 1));
}

TEST_F(FixedArraysTest, LowerTriangularMatrixIndexing) {
  EXPECT_EQ(10, (FixedLowerTriangularMatrix<double, 4>::size()));
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
}

TEST_F(FixedArraysTest, UpperTriangularMatrixIndexing) {
  EXPECT_EQ(10, (FixedUpperTriangularMatrix<double, 4>::size()));
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

}  // namespace numerics
}  // namespace principia
