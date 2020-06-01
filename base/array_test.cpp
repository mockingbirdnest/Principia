
#include "base/array.hpp"

#include <string>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::Not;
using ::testing::SizeIs;

namespace principia {
namespace base {

TEST(ArrayTest, BoundedArray0) {
  BoundedArray<double, 3> a{};
  EXPECT_THAT(a, AllOf(IsEmpty(), SizeIs(0)));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(0.0));
}

TEST(ArrayTest, BoundedArray1) {
  BoundedArray<double, 3> a{1.0};
  EXPECT_THAT(a, AllOf(Not(IsEmpty()), SizeIs(1)));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(1.0));
}

TEST(ArrayTest, BoundedArray2) {
  BoundedArray<double, 3> a{1.0, 10.0};
  EXPECT_THAT(a, AllOf(Not(IsEmpty()), SizeIs(2)));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(11.0));
}

TEST(ArrayTest, BoundedArray3) {
  BoundedArray<double, 3> a{1.0, 10.0, 100.0};
  EXPECT_THAT(a, AllOf(Not(IsEmpty()), SizeIs(3)));
  double sum = 0.0;
  for (double const value : a) {
    sum += value;
  }
  EXPECT_THAT(sum, Eq(111.0));
}

TEST(ArrayTest, BoundedArrayPushBack) {
  BoundedArray<double, 3> a{};
  BoundedArray<double, 5> const b{1.0, 2.0, 3.0};
  for (double const value : b) {
    a.push_back(value);
  }
  EXPECT_THAT(a, ElementsAre(1.0, 2.0, 3.0));
}

TEST(ArrayTest, BoundedArrayReverse) {
  BoundedArray<double, 3> a{};
  BoundedArray<double, 5> const b{1.0, 2.0, 3.0};
  for (auto it = b.rbegin(); it != b.rend(); ++it) {
    a.push_back(*it);
  }
  EXPECT_THAT(a, ElementsAre(3.0, 2.0, 1.0));
}

TEST(ArrayTest, Return) {
  std::string s;
  auto fn = [s]() -> BoundedArray<std::string, 3> {
    return {s};
  };
}

// This test only compiles if the constructor correctly uses |std::forward|.
TEST(ArrayTest, Move) {
  [[maybe_unused]] auto fn =
      [](std::unique_ptr<int> x,
         std::unique_ptr<int> y) -> BoundedArray<std::unique_ptr<int>, 3> {
    return {std::move(x), std::move(y)};
  };
}

}  // namespace base
}  // namespace principia
