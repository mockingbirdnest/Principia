#include "base/not_null.hpp"

#include <memory>
#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

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
    not_null<int *> int_ptr = check_not_null(null_int_ptr);
  }, "Check failed: .* != nullptr");
  EXPECT_DEATH({
    std::unique_ptr<int> null_int_ptr;
    not_null<std::unique_ptr<int>> int_ptr =
        check_not_null(std::move(null_int_ptr));
  }, "Check failed: .* != nullptr");
}

TEST_F(NotNullDeathTest, UniqueAssignment) {
  not_null<std::unique_ptr<int>> int_ptr1 =
      check_not_null(std::make_unique<int>(3));
  not_null<std::unique_ptr<int>> int_ptr2 = std::move(int_ptr1);
}

TEST_F(NotNullDeathTest, Assignment) {
  std::unique_ptr<int> int_ptr = std::make_unique<int>(3);
  not_null<int*> int_ptr1 = check_not_null(int_ptr.get());
  not_null<int*> int_ptr2 = int_ptr1;
}


}  // namespace base
}  // namespace principia
