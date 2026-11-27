import std;

import "absl/container/btree_map.h";
import "absl/container/btree_set.h";
import "absl/container/fixed_array.h";
import "absl/container/flat_hash_map.h";
import "absl/container/flat_hash_set.h";
import "absl/container/inlined_vector.h";
import "absl/container/node_hash_map.h";
import "absl/container/node_hash_set.h";
import "gmock/gmock.h";
import "gtest/gtest.h";

import principia.base.malloc_allocator;

// "macos_allocator_replacement.hpp" should be automagically included.

#ifdef OS_MACOSX

namespace principia {
namespace base {

using namespace principia::base::_malloc_allocator;

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
