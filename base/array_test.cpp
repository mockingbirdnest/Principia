
#include "base/array.hpp"

#include <string>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using ::testing::Eq;

namespace principia {
namespace base {

TEST(ArrayTest, BoundedArray0) {
  BoundedArray<double, 3> a{};
  EXPECT_THAT(a.size(), Eq(0));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(0.0));
}

TEST(ArrayTest, BoundedArray1) {
  BoundedArray<double, 3> a{1.0};
  EXPECT_THAT(a.size(), Eq(1));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(1.0));
}

TEST(ArrayTest, BoundedArray2) {
  BoundedArray<double, 3> a{1.0, 10.0};
  EXPECT_THAT(a.size(), Eq(2));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(11.0));
}

TEST(ArrayTest, BoundedArray3) {
  BoundedArray<double, 3> a{1.0, 10.0, 100.0};
  EXPECT_THAT(a.size(), Eq(3));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(111.0));
}

TEST(ArrayTest, Return) {
  std::string s;
  auto fn = [s]() -> BoundedArray<std::string, 3> {
    return {s};
  };
}

// This test only compiles if the constructor correctly uses |std::forward|.
TEST(ArrayTest, Move) {
  auto fn =
      [](std::unique_ptr<int> x,
         std::unique_ptr<int> y) -> BoundedArray<std::unique_ptr<int>, 3> {
    return {std::move(x), std::move(y)};
  };
}

}  // namespace base
}  // namespace principia
