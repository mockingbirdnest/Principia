#include "base/allocated_by.hpp"

#include <memory>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using testing::IsNull;
using testing::Not;
using testing::Pointee;

TEST(AllocatedByTest, NullablePointer) {
  // Test that AllocatedBy satisfies the NullablePointer C++ named
  // requirement.

  AllocatedBy<std::allocator<int>> a(nullptr);
  EXPECT_THAT(a, IsNull());

  AllocatedBy<std::allocator<int>> b = nullptr;
  EXPECT_THAT(b, IsNull());

  EXPECT_THAT(AllocatedBy<std::allocator<int>>(nullptr), IsNull());

  std::vector<int> v(1);
  auto c = AssertAllocatedBy<std::allocator<int>>(&v[0]);
  EXPECT_THAT(c, Not(IsNull()));
  EXPECT_THAT(c = nullptr, IsNull());
  EXPECT_THAT(c, IsNull());

  auto p = AssertAllocatedBy<std::allocator<int>>(&v[0]);
  AllocatedBy<std::allocator<int>> q = nullptr;
  EXPECT_FALSE(p == q);
  EXPECT_TRUE(p != q);
  EXPECT_FALSE(p == nullptr);
  EXPECT_FALSE(nullptr == p);
  EXPECT_TRUE(q == nullptr);
  EXPECT_TRUE(nullptr == q);
  EXPECT_TRUE(p != nullptr);
  EXPECT_TRUE(nullptr != p);
  EXPECT_FALSE(q != nullptr);
  EXPECT_FALSE(nullptr != q);
}

TEST(AllocatedByTest, DefaultConstructor) {
  EXPECT_THAT(AllocatedBy<std::allocator<int>>(), IsNull());
}

TEST(AllocatedByTest, Dereference) {
  std::vector<std::string> v = {"foo"};
  AllocatedBy<std::allocator<std::string>> p =
      AssertAllocatedBy<std::allocator<std::string>>(&v[0]);

  EXPECT_EQ(*p, "foo");
  EXPECT_EQ(p->size(), 3);
}

TEST(AllocatedByTest, Lifecycle) {
  std::vector<int> v(1);
  int* ptr = &v[0];
  auto asserted_ptr = AssertAllocatedBy<std::allocator<int>>(ptr);
  EXPECT_EQ(asserted_ptr, ptr);
}

class Base {};
class Sub : public Base {};

TEST(AllocatedByTest, Upcasting) {
  // Should compile.
  AllocatedBy<std::allocator<Sub>> a = nullptr;
  AllocatedBy<std::allocator<Base>> b = a;
  EXPECT_THAT(b, IsNull());
}

}  // namespace base

}  // namespace principia
