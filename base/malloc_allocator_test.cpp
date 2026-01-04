#include "base/malloc_allocator.hpp"

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::IsNull;
using ::testing::Not;
using namespace principia::base::_malloc_allocator;

struct alignas(32) S {};

TEST(MallocAllocatorTest, Allocate) {
  int* p;
  p = MallocAllocator<int>().allocate(1);
  EXPECT_THAT(p, Not(IsNull()));
  *p = 123;
#if PRINCIPIA_COMPILER_MSVC
    _aligned_free(p);
#else
    free(p);
#endif

  p = MallocAllocator<int>().allocate(10);
  EXPECT_THAT(p, Not(IsNull()));
  p[9] = 123;
#if PRINCIPIA_COMPILER_MSVC
    _aligned_free(p);
#else
    free(p);
#endif
}

TEST(MallocAllocatorTest, Deallocate) {
#if PRINCIPIA_COMPILER_MSVC
    int* p = static_cast<int*>(_aligned_malloc(sizeof(int), alignof(int)));
#else
    int* p = static_cast<int*>(malloc(sizeof(int)));
#endif
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

TEST(MallocAllocatorTest, Alignment) {
  S* p;
  for (int i = 0; i < 10; ++i) {
    p = MallocAllocator<S>().allocate(1);
    EXPECT_EQ(0, reinterpret_cast<intptr_t>(p) % alignof(S));
  }
}

}  // namespace base
}  // namespace principia
