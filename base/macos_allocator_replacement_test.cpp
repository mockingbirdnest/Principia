
// "macos_allocator_replacement.hpp" should be automagically included.

#include <deque>
#include <list>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "absl/container/fixed_array.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/container/node_hash_map.h"
#include "absl/container/node_hash_set.h"
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

  // Abseil
  static_assert(AllocatorIs<absl::btree_set<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<absl::btree_map<int, int>,
                            MallocAllocator<std::pair<const int, int>>>());
  static_assert(AllocatorIs<absl::flat_hash_set<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<absl::flat_hash_map<int, int>,
                            MallocAllocator<std::pair<const int, int>>>());
  static_assert(AllocatorIs<absl::node_hash_set<int>, MallocAllocator<int>>());
  static_assert(AllocatorIs<absl::node_hash_map<int, int>,
                            MallocAllocator<std::pair<const int, int>>>());
  static_assert(
      AllocatorIs<absl::FixedArray<int, 1024>, MallocAllocator<int>>());
  static_assert(
      AllocatorIs<absl::InlinedVector<int, 4>, MallocAllocator<int>>());
}

}  // namespace base
}  // namespace principia

#endif
