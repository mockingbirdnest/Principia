
#include "base/array.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace base {

TEST(ArrayTest, BoundedArray0) {
  BoundedArray<double, 3> a{{}, 0};
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_EQ(0.0, sum);
}

TEST(ArrayTest, BoundedArray1) {
  BoundedArray<double, 3> a{{1.0}, 1};
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_EQ(1.0, sum);
}

TEST(ArrayTest, BoundedArray2) {
  BoundedArray<double, 3> a{{1.0, 10.0}, 2};
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_EQ(11.0, sum);
}

TEST(ArrayTest, BoundedArray3) {
  BoundedArray<double, 3> a{{1.0, 10.0, 100.0}, 3};
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_EQ(111.0, sum);
}

}  // namespace base
}  // namespace principia
