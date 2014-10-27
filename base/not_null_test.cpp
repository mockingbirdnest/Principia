#include "base/not_null.hpp"

#include <memory>
#include <utility>
#include <string>

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

TEST_F(NotNullTest, Move) {
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

TEST_F(NotNullTest, Copy) {
  std::unique_ptr<int> const owner_of_three = std::make_unique<int>(3);
  std::unique_ptr<int> const owner_of_five = std::make_unique<int>(5);
  not_null<int*> const int_ptr1 = check_not_null(owner_of_three.get());
  not_null<int*> int_ptr2 = int_ptr1;
  not_null<int*> const int_ptr3(int_ptr2);
  EXPECT_THAT(int_ptr1, Eq(owner_of_three.get()));
  EXPECT_THAT(int_ptr2, Eq(owner_of_three.get()));
  EXPECT_THAT(int_ptr3, Eq(owner_of_three.get()));
  EXPECT_THAT(*int_ptr1, Eq(3));
  EXPECT_THAT(*int_ptr2, Eq(3));
  EXPECT_THAT(*int_ptr3, Eq(3));
  int_ptr2 = check_not_null(owner_of_five.get());
  EXPECT_THAT(*int_ptr2, Eq(5));
  int_ptr2 = std::move(int_ptr3);
  EXPECT_THAT(*int_ptr2, Eq(3));
  EXPECT_THAT(*int_ptr3, Eq(3));
}

TEST_F(NotNullTest, CheckNotNull) {
  std::unique_ptr<int> owner_int = std::make_unique<int>(3);
  int* const constant_access_int = owner_int.get();
  not_null<int*> const constant_not_null_access_int =
      check_not_null(constant_access_int);
  check_not_null(constant_not_null_access_int);
  not_null<std::unique_ptr<int>> not_null_owner_int =
      check_not_null(std::move(owner_int));
  check_not_null(std::move(not_null_owner_int));
}

TEST_F(NotNullTest, Booleans) {
  not_null<std::unique_ptr<int>> const pointer =
      check_not_null(std::make_unique<int>(3));
  EXPECT_TRUE(pointer);
  EXPECT_TRUE(pointer != nullptr);
  EXPECT_FALSE(pointer == nullptr);
}

TEST_F(NotNullTest, ImplicitConversions) {
  not_null<std::unique_ptr<int>> not_null_owner_int =
      check_not_null(std::make_unique<int>(3));
  not_null<int const*> const constant_not_null_access_int =
      not_null_owner_int.get();
  // Copy constructor.
  not_null<int const*> not_null_access_constant_int =
      constant_not_null_access_int;
  // Copy assignment.
  not_null_access_constant_int = not_null_owner_int.get();
  // Move constructor.
  not_null<std::unique_ptr<int const>> not_null_owner_constant_int =
      check_not_null(std::make_unique<int>(5));
  // Move assigment.
  not_null_owner_constant_int = std::move(not_null_owner_int);
}

TEST_F(NotNullTest, Arrow) {
  not_null<std::unique_ptr<std::string>> not_null_owner_string =
      check_not_null(std::make_unique<std::string>("-"));
  not_null_owner_string->append(">");
  EXPECT_THAT(*not_null_owner_string, Eq("->"));
  not_null<std::string*> not_null_access_string = not_null_owner_string.get();
  not_null_access_string->insert(0, "operator");
  EXPECT_THAT(*not_null_access_string, Eq("operator->"));
}

}  // namespace base
}  // namespace principia
