#include "base/malloc_allocator.hpp"

#include <cstdlib>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::IsNull;
using ::testing::Not;
using namespace principia::base::_malloc_allocator;

TEST(MallocAllocatorTest, Allocate) {
  int* p;
  p = MallocAllocator<int>().allocate(1);
  EXPECT_THAT(p, Not(IsNull()));
  *p = 123;
  free(p);

  p = MallocAllocator<int>().allocate(10);
  EXPECT_THAT(p, Not(IsNull()));
  p[9] = 123;
  free(p);
}

TEST(MallocAllocatorTest, Deallocate) {
  int* p = static_cast<int*>(malloc(sizeof(int)));
  MallocAllocator<int>().deallocate(p, 1);
}

TEST(MallocAllocatorTest, RoundTrip) {
  int* p;
  p = MallocAllocator<int>().allocate(1);
  EXPECT_THAT(p, Not(IsNull()));
  *p = 123;
  MallocAllocator<int>().deallocate(p, 1);
}

TEST(MallocAllocatorTest, Conversion) {
  // This is required by some containers.
  MallocAllocator<int> foo;
  MallocAllocator<float> bar(foo);
}

}  // namespace base
}  // namespace principia
