#include "numerics/arrays.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace numerics {

class ArraysTest : public ::testing::Test {
 protected:
  ArraysTest()
    : m34_({-8, -6, -4, -7,
            -4, -10, 9, -5,
             6, -3, -2, -9}),
      v3_({10, 31, -47}),
      v4_({-3, -3, 1, 4}) {}

  FixedMatrix<double, 3, 4> m34_;
  FixedVector<double, 3> v3_;
  FixedVector<double, 4> v4_;
};

TEST_F(ArraysTest, Assignment) {
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

TEST_F(ArraysTest, Multiplication) {
  EXPECT_EQ(v3_, m34_ * v4_);
}

TEST_F(ArraysTest, Indexing) {
  EXPECT_EQ(31, v3_[1]);
  v3_[2] = -666;
  EXPECT_EQ(-666, v3_[2]);
}

}  // namespace numerics
}  // namespace principia
