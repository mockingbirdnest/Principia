#include "numerics/unbounded_arrays.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/transposed_view.hpp"
#include "quantities/elementary_functions.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::testing_utilities::_almost_equals;
using ::testing::Pointer;

class UnboundedArraysTest : public ::testing::Test {
 protected:
  UnboundedArraysTest()
    : u3_({6, -1, 12}),
      v3_({10, 31, -47}),
      v4_({-3, -3, 1, 4}),
      m34_(3, 4,
           {-8,  -6, -4, -7,
            -4, -10,  9, -5,
             6,  -3, -2, -9}),
      m23_(2, 3,
           {1, -2,  0,
            2,  3,  7}),
      n23_(2, 3,
           { 5, -1,  3,
            12, 13, -4}),
      l4_({ 1,
            2,  3,
            5,  8,  13,
            21, 34, 55, 89}),
      su4_({    1, 2,  3,
                   5,  8,
                      13
                        }),
      u4_({1, 2,  3,  5,
              8, 13, 21,
                 34, 55,
                     89}) {}

  UnboundedVector<double> u3_;
  UnboundedVector<double> v3_;
  UnboundedVector<double> v4_;
  UnboundedMatrix<double> m34_;
  UnboundedMatrix<double> m23_;
  UnboundedMatrix<double> n23_;
  UnboundedLowerTriangularMatrix<double> l4_;
  UnboundedStrictlyUpperTriangularMatrix<double> su4_;
  UnboundedUpperTriangularMatrix<double> u4_;
};

TEST_F(UnboundedArraysTest, Assignment) {
  {
    UnboundedVector<double> u2({1, 2});
    UnboundedVector<double> v2 = {{1, 2}};
    UnboundedVector<double> w2(2);
    w2 = {1, 2};
    EXPECT_EQ(u2, v2);
    EXPECT_EQ(u2, w2);
  }
  {
    UnboundedMatrix<double> l23(2, 3,
                                {1, 2, 3,
                                 4, 5, 6 });
    UnboundedMatrix<double> m23 = {2, 3,
                                   {1, 2, 3,
                                    4, 5, 6}};
    UnboundedMatrix<double> n23 = {2, 3,
                                   {0, 0, 0,
                                    0, 0, 0}};
    n23 = {1, 2, 3,
           4, 5, 6};
    EXPECT_EQ(l23, m23);
    EXPECT_EQ(l23, n23);
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
    n3 = {1,
          2, 3,
          4, 5, 6};
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
    n3 = {1, 2, 3,
             4, 5,
                6};
    EXPECT_EQ(l3, m3);
    EXPECT_EQ(l3, n3);
  }
}

TEST_F(UnboundedArraysTest, Norm) {
  EXPECT_EQ(35, TransposedView{v4_} * v4_);  // NOLINT
  EXPECT_EQ(Sqrt(35.0), v4_.Norm());
  EXPECT_EQ(35, v4_.Norm²());
  EXPECT_EQ(Sqrt(517.0), m34_.FrobeniusNorm());
}

TEST_F(UnboundedArraysTest, AdditiveGroups) {
  EXPECT_EQ((UnboundedVector<double>({-10, -31, 47})), -v3_);
  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     {-1,  2,  0,
                                      -2, -3, -7})), -m23_);

  EXPECT_EQ((UnboundedVector<double>({16, 30, -35})), u3_ + v3_);
  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     { 6, -3, 3,
                                      14, 16, 3})), m23_ + n23_);

  EXPECT_EQ((UnboundedVector<double>({-4, -32, 59})), u3_ - v3_);
  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     { -4,  -1, -3,
                                      -10, -10, 11})), m23_ - n23_);
}

TEST_F(UnboundedArraysTest, VectorSpaces) {
  EXPECT_EQ((UnboundedVector<double>({12, -2, 24})), 2 * u3_);
  EXPECT_EQ((UnboundedVector<double>({-30, -93, 141})), v3_ * -3);

  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     {2, -4,  0,
                                      4,  6, 14})), 2 * m23_);
  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     {-15,   3, -9,
                                      -36, -39, 12})), n23_ * -3);

  EXPECT_EQ((UnboundedVector<double>({2.5, 7.75, -11.75})), v3_ / 4);
  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     {-2.5,  0.5, -1.5,
                                        -6, -6.5,  2})), n23_ / -2);
}

TEST_F(UnboundedArraysTest, Algebra) {
  EXPECT_EQ(-535, TransposedView{u3_} * v3_);  // NOLINT
  EXPECT_EQ((UnboundedMatrix<double>(3, 4,
                                     { -30, -30,  10,   40,
                                       -93, -93,  31,  124,
                                       141, 141, -47, -188})),
             v3_ * TransposedView{v4_});
  EXPECT_EQ((UnboundedMatrix<double>(2, 4,
                                     { 0,  14, -22,   3,
                                      14, -63,   5, -92})),
            m23_ * m34_);
  auto m = UnboundedMatrix<double>({1, 2,
                                    3, 4});
  EXPECT_THAT(&(m *= UnboundedMatrix<double>({5, 6,
                                              7, 8})),
              Pointer(&m));
  EXPECT_EQ((UnboundedMatrix<double>({19, 22,
                                      43, 50})),
            m);
  EXPECT_EQ(v3_, m34_ * v4_);
  EXPECT_EQ((UnboundedVector<double>({-486, -229, 333, 198})),
            TransposedView{m34_} * v3_);  // NOLINT
  UnboundedMatrix<double> m43(TransposedView{m34_});
  EXPECT_EQ((UnboundedVector<double>({-486, -229, 333, 198})),
            m43 * v3_);
}

TEST_F(UnboundedArraysTest, VectorIndexing) {
  EXPECT_EQ(31, v3_[1]);
  v3_[2] = -666;
  EXPECT_EQ(-666, v3_[2]);
}

TEST_F(UnboundedArraysTest, MatrixIndexing) {
  EXPECT_EQ(9, m34_(1, 2));
  m34_(2, 1) = -666;
  EXPECT_EQ(-666, m34_(2, 1));
}

TEST_F(UnboundedArraysTest, LowerTriangularMatrixIndexing) {
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

  UnboundedLowerTriangularMatrix<double> const l4 = l4_;
  EXPECT_EQ(1, l4(0, 0));
}

TEST_F(UnboundedArraysTest, StrictlyUpperTriangularMatrixIndexing) {
  EXPECT_EQ(1, su4_(0, 1));
  EXPECT_EQ(2, su4_(0, 2));
  EXPECT_EQ(3, su4_(0, 3));
  EXPECT_EQ(5, su4_(1, 2));
  EXPECT_EQ(8, su4_(1, 3));
  EXPECT_EQ(13, su4_(2, 3));
  su4_(1, 3) = -666;
  EXPECT_EQ(-666, su4_(1, 3));

  UnboundedStrictlyUpperTriangularMatrix<double> const su4 = su4_;
  EXPECT_EQ(1, su4(0, 1));
}

TEST_F(UnboundedArraysTest, UpperTriangularMatrixIndexing) {
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

  UnboundedUpperTriangularMatrix<double> const u4 = u4_;
  EXPECT_EQ(1, u4(0, 0));
}

TEST_F(UnboundedArraysTest, Conversions) {
  EXPECT_EQ(UnboundedMatrix<double>(4, 4,
                                    {1,   0,  0,  0,
                                     2,   3,  0,  0,
                                     5,   8, 13,  0,
                                     21, 34, 55, 89}),
            UnboundedMatrix<double>(l4_));
  EXPECT_EQ(UnboundedMatrix<double>(4, 4,
                                    {0, 1, 2,  3,
                                     0, 0, 5,  8,
                                     0, 0, 0, 13,
                                     0, 0, 0,  0}),
            UnboundedMatrix<double>(su4_));
  EXPECT_EQ(UnboundedMatrix<double>(4, 4,
                                    {1, 2,  3,  5,
                                     0, 8, 13, 21,
                                     0, 0, 34, 55,
                                     0, 0,  0, 89}),
            UnboundedMatrix<double>(u4_));
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

TEST_F(UnboundedArraysTest, Transpose) {
  EXPECT_EQ(
      UnboundedMatrix<double>(4, 3,
                              {-8,  -4,  6,
                               -6, -10, -3,
                               -4,   9, -2,
                               -7,  -5, -9}),
      UnboundedMatrix<double>(TransposedView{m34_}));
  EXPECT_EQ(
      UnboundedUpperTriangularMatrix<double>({1, 2,  5, 21,
                                                 3,  8, 34,
                                                    13, 55,
                                                        89}),
      UnboundedUpperTriangularMatrix<double>(TransposedView{l4_}));
  EXPECT_EQ(
      UnboundedLowerTriangularMatrix<double>({1,
                                              2,  8,
                                              3, 13, 34,
                                              5, 21, 55, 89}),
      UnboundedLowerTriangularMatrix<double>(TransposedView{u4_}));
}

}  // namespace numerics
}  // namespace principia
