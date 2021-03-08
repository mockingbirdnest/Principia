#include "base/allocator_new.hpp"

#include <memory>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using testing::IsNull;
using testing::Not;
using testing::Pointee;

TEST(AllocatorNewTest, New) {
  int* p = new (AllocateWith<std::allocator<int>>{}) int(2);
  EXPECT_THAT(p, Pointee(2));
  std::allocator<int>().deallocate(p, 1);
}

TEST(AllocatorNewTest, Const) {
  int const* p = new (AllocateWith<std::allocator<int const>>{}) int const(2);
  EXPECT_THAT(p, Pointee(2));
  std::allocator<int const>().deallocate(p, 1);
}

TEST(AllocatorNewDeathTest, SizeMismatch) {
  EXPECT_DEATH(new (AllocateWith<std::allocator<uint8_t>>{}) uint16_t,
               "Number of bytes allocated must equal size of object.");
}

TEST(AllocatorNewTest, NonAllocatingPlacementNew) {
  // Check that we didn't somehow break non-allocating placement new.
  int a = 2;
  new (&a) int(3);
  EXPECT_EQ(a, 3);
}

}  // namespace base
}  // namespace principia
