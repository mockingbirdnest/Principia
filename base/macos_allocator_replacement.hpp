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

#include "base/macros.hpp"
#include "base/malloc_allocator.hpp"

#if !OS_MACOSX
#error Only include this file for macOS
#endif

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

}  // namespace std
}  // namespace principia
