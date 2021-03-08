
#pragma once

#include <cstddef>
#include <cstdlib>
#include <type_traits>

namespace principia {
namespace base {

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

  T* allocate(size_t n) {
    return static_cast<T*>(calloc(n, sizeof(T)));
  }

  void deallocate(T* p, size_t n) {
    free(const_cast<std::remove_const_t<T>*>(p));
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

}  // namespace base
}  // namespace principia
