
#pragma once

#include <memory>
#include <type_traits>

#include "glog/logging.h"

namespace principia {
namespace base {
// Tag type for passing to allocator new.
template<typename Alloc>
struct AllocateWith {
  static_assert(std::is_empty<Alloc>::value, "|Alloc| must be stateless");
};
}  // namespace base
}  // namespace principia

// Placement new operator using allocators for allocation.
// Should only be used for allocating single objects (for arrays and such, use
// the allocator directly).
template<typename Alloc>
void* operator new(std::size_t count,
                   principia::base::AllocateWith<Alloc>) noexcept {
  constexpr size_t size =
      sizeof(typename std::allocator_traits<Alloc>::value_type);
  CHECK_EQ(count, size)
      << "Number of bytes allocated must equal size of object.";
  Alloc alloc;
  auto p = std::allocator_traits<Alloc>::allocate(alloc, 1);
  using non_const_T =
      std::remove_const_t<typename std::allocator_traits<Alloc>::value_type>;
  return const_cast<non_const_T*>(p);
}

// The corresponding delete operator.
// Note that this is only called when a class constructor throws an exception.
template<typename Alloc>
void operator delete(void* ptr, principia::base::AllocateWith<Alloc>) noexcept {
  using pointer = typename std::allocator_traits<Alloc>::pointer;
  Alloc alloc;
  std::allocator_traits<Alloc>::deallocate(alloc, static_cast<pointer>(ptr), 1);
}
