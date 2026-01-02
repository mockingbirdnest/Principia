#pragma once

#include <cstddef>
#include <cstdlib>

#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.

namespace principia {
namespace base {
namespace _malloc_allocator {
namespace internal {

// An allocator (for use with containers such as `std::vector`) that uses malloc
// and free for memory management instead of global new and delete. The purpose
// of this allocator is to enable use of the system allocator even when global
// new and delete have been overridden.
template<typename T>
class MallocAllocator {
 public:
  using value_type = T;
  using pointer = T*;
  using const_pointer = T const*;
  using reference = T&;
  using const_reference = T const&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  // Constructors. This is a stateless struct so these don't do anything (but
  // they are required by some containers).
  MallocAllocator() {}
  MallocAllocator(const MallocAllocator& other) = default;
  template<typename U>
  MallocAllocator(const MallocAllocator<U>& other) {}

  T* allocate(std::size_t n) {
#if PRINCIPIA_COMPILER_MSVC
    return static_cast<T*>(_aligned_malloc(n * sizeof(T), alignof(T)));
#else
    // On macOS, aligned_alloc returns null for small alignments; use malloc if
    // that is aligned enough.
    if constexpr (alignof(T) > alignof(std::max_align_t)) {
      return static_cast<T*>(std::aligned_alloc(alignof(T), n * sizeof(T)));
    } else {
      return static_cast<T*>(std::malloc(n * sizeof(T)));
    }
#endif
  }

  void deallocate(T* p, std::size_t n) {
#if PRINCIPIA_COMPILER_MSVC
    _aligned_free(p);
#else
    free(p);
#endif
  }
};

// MallocAllocators are equal regardless of type.
template<typename T1, typename T2>
constexpr bool operator==(const MallocAllocator<T1>&,
                          const MallocAllocator<T2>&) {
  return true;
}

template<typename T1, typename T2>
constexpr bool operator!=(const MallocAllocator<T1>&,
                          const MallocAllocator<T2>&) {
  return false;
}

}  // namespace internal

using internal::MallocAllocator;
using internal::operator!=;
using internal::operator==;

}  // namespace _malloc_allocator
}  // namespace base
}  // namespace principia
