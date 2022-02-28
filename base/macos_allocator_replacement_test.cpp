
// "macos_allocator_replacement.hpp" should be automagically included.

#include <deque>
#include <list>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/macros.hpp"
#include "base/malloc_allocator.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#ifdef OS_MACOSX

namespace principia {
namespace base {

// AllocatorIs<container, alloc>() returns true iff container's allocator is
// alloc.
template <typename T, typename Allocator>
constexpr bool AllocatorIs() {
  return std::is_same<typename T::allocator_type, Allocator>::value;
}

// Test that the default allocators for various classes are overridden.
TEST(PrincipiaMallocAllocatorTest, DefaultAllocators) {
  // STL
  static_assert(AllocatorIs<std::vector<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<std::deque<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<std::list<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<std::set<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<std::map<int, int>,
                            MallocAllocator<std::pair<const int, int>>>());
  static_assert(AllocatorIs<std::multiset<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<std::multimap<int, int>,
                            MallocAllocator<std::pair<const int, int>>>());
}

}  // namespace base
}  // namespace principia

#endif
