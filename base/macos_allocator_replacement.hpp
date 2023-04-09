// This file contains aliases of some STL and Abseil containers (such as
// `std::vector`) in the `principia::std` and `principia::absl` namespaces with
// their default allocators overridden to base new one defined here. The result
// of this is that if you write `std::vector foo` inside the principia
// namespace, it will resolve to a vector with the allocator defined here.
//
// The intended purpose of this file is to be used on macOS builds. This is
// necessary because Unity provides overrides of `new` and `delete` which have
// terrible performance on macOS. It should be included via compiler argument.
// Attempting to include it on non-macOS builds will result in an error.

#pragma once

#include <deque>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <set>
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

#if !OS_MACOSX
#error Only include this file for macOS
#endif

namespace principia {
namespace std {

using namespace ::std;

template <typename T>
using allocator = ::principia::base::_malloc_allocator::MallocAllocator<T>;

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

namespace absl {

using namespace ::absl;

template <typename Key, typename Compare = ::std::less<Key>,
          typename Allocator = std::allocator<Key>>
using btree_set = ::absl::btree_set<Key, Compare, Allocator>;

template <typename Key, typename Value, typename Compare = ::std::less<Key>,
          typename Allocator = std::allocator<::std::pair<const Key, Value>>>
using btree_map = ::absl::btree_map<Key, Value, Compare, Allocator>;

template <typename T,
          typename Hash = ::absl::container_internal::hash_default_hash<T>,
          typename Eq = ::absl::container_internal::hash_default_eq<T>,
          typename Allocator = std::allocator<T>>
using flat_hash_set = ::absl::flat_hash_set<T, Hash, Eq, Allocator>;

template <typename K, typename V,
          typename Hash = ::absl::container_internal::hash_default_hash<K>,
          typename Eq = ::absl::container_internal::hash_default_eq<K>,
          typename Allocator = std::allocator<::std::pair<const K, V>>>
using flat_hash_map = ::absl::flat_hash_map<K, V, Hash, Eq, Allocator>;

template <typename T,
          typename Hash = ::absl::container_internal::hash_default_hash<T>,
          typename Eq = ::absl::container_internal::hash_default_eq<T>,
          typename Allocator = std::allocator<T>>
using node_hash_set = ::absl::node_hash_set<T, Hash, Eq, Allocator>;

template <typename K, typename V,
          typename Hash = ::absl::container_internal::hash_default_hash<K>,
          typename Eq = ::absl::container_internal::hash_default_eq<K>,
          typename Allocator = std::allocator<::std::pair<const K, V>>>
using node_hash_map = ::absl::node_hash_map<K, V, Hash, Eq, Allocator>;

template <typename T, size_t N = ::absl::kFixedArrayUseDefault,
          typename Allocator = std::allocator<T>>
using FixedArray = ::absl::FixedArray<T, N, Allocator>;

template <typename T, size_t N, typename Allocator = std::allocator<T>>
using InlinedVector = ::absl::InlinedVector<T, N, Allocator>;

}  // namespace absl
}  // namespace principia
