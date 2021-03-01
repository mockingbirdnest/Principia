// This file contains two things:
//
// 1. An allocator (for use with containers such as `std::vector`) that uses
//    malloc and free for memory management instead of global new and delete.
//    The purpose of this allocator is to enable use of the system allocator
//    even when global new and delete have been overridden.
//
// 2. Aliases of some STL containers (such as `std::vector`) in the
//    `principia::std` namespace with their default allocators overridden to the
//    new one defined here. The result of this is that if you write `std::vector
//    foo` inside the principia namespace, it will resolve to a vector with the
//    allocator defined here.
//
// That is, if you include this file and write code in the `principia`
// namespace, it will automatically use the system allocator even when new and
// delete have been overridden.
//
// The intended purpose of this file is to be used on macOS builds. This is
// necessary because Unity provides overrides of `new` and `delete` which have
// terrible performance on macOS.
//
// This file is placed in the "magic" directory because it is not intended to be
// included directly. Instead, it should be included "magically" via compiler
// argument.

#pragma once

#include <cstddef>
#include <cstdlib>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace principia {
namespace magic {

// Allocator which uses malloc/free instead of new/delete. This is useful when
// undesired global new and delete overloads are present.
template <typename T>
class MallocAllocator {
 public:
  using value_type = T;
  using pointer = T*;
  using const_pointer = T const*;
  using reference = T&;
  using const_reference = T const&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  T* allocate(size_t n) { return static_cast<T*>(calloc(n, sizeof(T))); }

  void deallocate(T* p, size_t n) { free(p); }
};

// MallocAllocators are equal regardless of type.
template <typename T1, typename T2>
constexpr bool operator==(const MallocAllocator<T1>&,
                          const MallocAllocator<T2>&) {
  return true;
}

template <typename T1, typename T2>
constexpr bool operator!=(const MallocAllocator<T1>&,
                          const MallocAllocator<T2>&) {
  return false;
}

}  // namespace magic

// We alias some things in namespace principia::std so that
// Principia code using std::vector and the like automatically use
// MallocAllocator.
namespace std {

using namespace ::std;

template <typename T>
using allocator = ::principia::magic::MallocAllocator<T>;

template <typename T, typename Allocator = allocator<T>>
using vector = ::std::vector<T, Allocator>;

template <typename T, typename Allocator = allocator<T>>
using deque = ::std::deque<T, Allocator>;

template <typename T, typename Allocator = allocator<T>>
using list = ::std::list<T, Allocator>;

template <typename Key, typename Compare = std::less<Key>,
          typename Allocator = allocator<Key>>
using set = ::std::set<Key, Compare, Allocator>;

template <typename Key, typename T, typename Compare = std::less<Key>,
          typename Allocator = allocator<std::pair<const Key, T>>>
using map = ::std::map<Key, T, Compare, Allocator>;

template <typename Key, typename Compare = std::less<Key>,
          typename Allocator = allocator<Key>>
using multiset = ::std::multiset<Key, Compare, Allocator>;

template <typename Key, typename T, typename Compare = std::less<Key>,
          typename Allocator = allocator<std::pair<const Key, T>>>
using multimap = ::std::multimap<Key, T, Compare, Allocator>;

}  // namespace std

}  // namespace principia
