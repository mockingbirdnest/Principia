
// "macos_allocator_replacement.hpp" should be automagically included.
#include "macos_allocator_replacement.hpp"

#include <deque>
#include <list>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/allocator_deleter.hpp"
#include "base/macros.hpp"
#include "base/malloc_allocator.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#ifdef OS_MACOSX

// namespace principia {
// namespace std {

// template<typename T,
//          typename Deleter =
//              ::std::conditional_t<!::std::is_array<T>::value,
//                                   base::AllocatorDeleter<allocator<T>>,
//                                   ::std::default_delete<T>>>
// using unique_ptr = ::std::unique_ptr<T, Deleter>;

// // C++ doesn't support function aliases, so we wrap make_unique instead.
// template<typename T, typename... Args>
// inline unique_ptr<T> make_unique(Args&&... args) {
//   return base::MakeUniqueWithAllocator<allocator<T>>(
//       std::forward<Args>(args)...);
// }

// }  // namespace std
// }  // namespace principia

namespace principia {
namespace base {

using testing::IsNull;
using testing::Not;
using testing::Pointee;

// AllocatorIs<container, alloc>() returns true iff container's allocator is
// alloc.
template<typename T, typename Allocator>
constexpr bool AllocatorIs() {
  return std::is_same<typename T::allocator_type, Allocator>::value;
}

// Test that the default allocators for various classes are overridden.
TEST(MacosAllocatorReplacementTest, DefaultAllocators) {
  // STL
  EXPECT_TRUE((AllocatorIs<std::vector<int>, MallocAllocator<int>>()));
  EXPECT_TRUE((AllocatorIs<std::deque<int>, MallocAllocator<int>>()));
  EXPECT_TRUE((AllocatorIs<std::list<int>, MallocAllocator<int>>()));
  EXPECT_TRUE((AllocatorIs<std::set<int>, MallocAllocator<int>>()));
  EXPECT_TRUE((AllocatorIs<std::map<int, int>,
                           MallocAllocator<std::pair<const int, int>>>()));
  EXPECT_TRUE((AllocatorIs<std::multiset<int>, MallocAllocator<int>>()));
  EXPECT_TRUE((AllocatorIs<std::multimap<int, int>,
                           MallocAllocator<std::pair<const int, int>>>()));
}

TEST(MacosAllocatorReplacementTest, UniquePtr) {
  // For non-array types, MallocAllocator should be used.
  EXPECT_TRUE((std::is_same<std::unique_ptr<int>::deleter_type,
                            AllocatorDeleter<MallocAllocator<int>>>::value));

  // For array types, ::std::default_deleter should be used.
  EXPECT_TRUE((std::is_same<std::unique_ptr<int[]>::deleter_type,
                            ::std::default_delete<int[]>>::value));

  // std::make_unique should also work appropriately.
  std::unique_ptr<int> p = std::make_unique<int>(2);
  EXPECT_THAT(p, Pointee(2));

  std::unique_ptr<int[]> arr = std::make_unique<int[]>(3);
  EXPECT_THAT(arr, Not(IsNull()));
}

}  // namespace base
}  // namespace principia

#endif
