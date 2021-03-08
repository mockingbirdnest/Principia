// This file contains aliases of some STL containers (such as `std::vector`) in
// the `principia::std` namespace with their default allocators overridden to
// base new one defined here. The result of this is that if you write
// `std::vector foo` inside the principia namespace, it will resolve to a vector
// with the allocator defined here.
//
// The intended purpose of this file is to be used on macOS builds. This is
// necessary because Unity provides overrides of `new` and `delete` which have
// terrible performance on macOS. It should be included via compiler argument.
// Attempting to include it on non-macOS builds will result in an error.

#pragma once

#include <deque>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "base/allocator_deleter.hpp"
#include "base/macros.hpp"
#include "base/malloc_allocator.hpp"

#if !OS_MACOSX
#error Only include this file for macOS
#endif

// Forward declare google::compression::Compressor.
namespace google {
namespace compression {
class Compressor;
}  // namespace compression
}  // namespace google

namespace principia {
namespace std {

using namespace ::std;

template<typename T>
using allocator = ::principia::base::MallocAllocator<T>;

template<typename T, typename Allocator = allocator<T>>
using vector = ::std::vector<T, Allocator>;

template<typename T, typename Allocator = allocator<T>>
using deque = ::std::deque<T, Allocator>;

template<typename T, typename Allocator = allocator<T>>
using list = ::std::list<T, Allocator>;

template<typename Key,
         typename Compare = std::less<Key>,
         typename Allocator = allocator<Key>>
using set = ::std::set<Key, Compare, Allocator>;

template<typename Key,
         typename T,
         typename Compare = std::less<Key>,
         typename Allocator = allocator<std::pair<const Key, T>>>
using map = ::std::map<Key, T, Compare, Allocator>;

template<typename Key,
         typename Compare = std::less<Key>,
         typename Allocator = allocator<Key>>
using multiset = ::std::multiset<Key, Compare, Allocator>;

template<typename Key,
         typename T,
         typename Compare = std::less<Key>,
         typename Allocator = allocator<std::pair<const Key, T>>>
using multimap = ::std::multimap<Key, T, Compare, Allocator>;

// unique_ptr is overridden for all classes except arrays and
// google::compression::Compressor.
template<typename T,
         typename Deleter = ::std::conditional_t<
             (::std::is_array<T>::value ||
              ::std::is_same<T, ::google::compression::Compressor>::value),
             ::std::default_delete<T>,
             base::AllocatorDeleter<allocator<T>>>>
using unique_ptr = ::std::unique_ptr<T, Deleter>;

// C++ doesn't support function aliases, so we wrap make_unique instead.

// Non-array variant of make_unique.
template<typename T, typename... Args>
inline unique_ptr<std::enable_if_t<!std::is_array<T>::value, T>> make_unique(
    Args&&... args) {
  return base::MakeUniqueWithAllocator<allocator<T>>(
      std::forward<Args>(args)...);
}

// Array variant of make_unique.
template<typename T>
inline unique_ptr<std::enable_if_t<std::is_array<T>::value, T>> make_unique(
    std::size_t size) {
  return ::std::make_unique<T>(size);
}

}  // namespace std
}  // namespace principia
