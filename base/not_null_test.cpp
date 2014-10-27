#include "base/not_null.hpp"

#include <memory>
#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using testing::Eq;

namespace principia {
namespace base {

class NotNullTest : public testing::Test {};

using NotNullDeathTest = NotNullTest;

TEST_F(NotNullDeathTest, DeathByNullptr) {
  EXPECT_DEATH({
    int* const null_int_ptr = nullptr;
    not_null<int*> int_ptr(null_int_ptr);
  }, "Check failed: .* != nullptr");

  EXPECT_DEATH({
    std::unique_ptr<int> null_int_ptr;
    not_null<std::unique_ptr<int>> int_ptr(std::move(null_int_ptr));
  }, "Check failed: .* != nullptr");

  EXPECT_DEATH({
    int* const null_int_ptr = nullptr;
    check_not_null(null_int_ptr);
  }, "Check failed: .* != nullptr");

  EXPECT_DEATH({
    std::unique_ptr<int> null_int_ptr;
    check_not_null(std::move(null_int_ptr));
  }, "Check failed: .* != nullptr");

}

TEST_F(NotNullDeathTest, Move) {
  not_null<std::unique_ptr<int>> int_ptr1 =
      check_not_null(std::make_unique<int>(3));
  EXPECT_THAT(*(std::unique_ptr<int> const&)int_ptr1, Eq(3));
  not_null<std::unique_ptr<int>> int_ptr2 = std::move(int_ptr1);
  EXPECT_THAT(*int_ptr2, Eq(3));
  not_null<std::unique_ptr<int>> int_ptr3(std::move(int_ptr2));
  EXPECT_THAT(*int_ptr3, Eq(3));
  int_ptr2 = check_not_null(std::make_unique<int>(5));
  EXPECT_THAT(*int_ptr2, Eq(5));
  int_ptr2 = std::move(int_ptr3);
  EXPECT_THAT(*int_ptr2, Eq(3));
  // This checks an implementation detail, namely that move assignment is
  // implemented as a swap.
  EXPECT_THAT(*int_ptr3, Eq(5));
}

TEST_F(NotNullDeathTest, Copy) {
  std::unique_ptr<int> owner_of_three = std::make_unique<int>(3);
  std::unique_ptr<int> owner_of_five = std::make_unique<int>(5);
  not_null<int*> int_ptr1 = check_not_null(owner_of_three.get());
  not_null<int*> int_ptr2 = int_ptr1;
  not_null<int*> int_ptr3(int_ptr2);
  EXPECT_THAT(*int_ptr1, Eq(3));
  EXPECT_THAT(*int_ptr2, Eq(3));
  EXPECT_THAT(*int_ptr3, Eq(3));
  int_ptr2 = check_not_null(owner_of_five.get());
  EXPECT_THAT(*int_ptr2, Eq(5));
  int_ptr2 = std::move(int_ptr3);
  EXPECT_THAT(*int_ptr2, Eq(3));
  EXPECT_THAT(*int_ptr3, Eq(3));
}


}  // namespace base
}  // namespace principia
