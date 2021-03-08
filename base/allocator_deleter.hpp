#pragma once

#include <memory>
#include <type_traits>

#include "base/allocated_by.hpp"
#include "base/allocator_new.hpp"

namespace principia {
namespace base {

// Deleter suitible for passing to std::unique_ptr that deletes using an
// allocator instead of ::delete.
template<typename Alloc>
class AllocatorDeleter {
 public:
  using traits = std::allocator_traits<Alloc>;
  using pointer = AllocatedBy<Alloc>;
  using value_type = typename traits::value_type;

  // We are only interested in stateless allocators.
  static_assert(std::is_empty<Alloc>::value);

  // Allocator deallocate requires the array size to be known so we can't handle
  // array types.
  static_assert(!std::is_array<typename traits::value_type>::value,
                "Array types are not supported.");

  constexpr AllocatorDeleter() = default;

  // Implicitly convertible from allocators of compatible types.
  // TODO(rnlahaye): check that allocator is compatible.
  template<typename B,
           typename = typename std::enable_if_t<
               std::is_convertible<typename std::allocator_traits<B>::pointer,
                                   typename traits::pointer>::value>>
  AllocatorDeleter(const AllocatorDeleter<B>& d) {}

  // Delete the managed object.
  void operator()(value_type* ptr) {
    // std::remove_const_t<value_type>* non_const_ptr =
    //     const_cast<std::remove_const_t<value_type>*>(ptr);
    auto non_const_ptr = ptr;
    Alloc alloc;
    traits::destroy(alloc, typename traits::pointer(non_const_ptr));
    traits::deallocate(alloc, non_const_ptr, 1);
  }
};

template<typename Alloc,
         typename T = typename std::allocator_traits<Alloc>::value_type,
         typename... Args>
std::unique_ptr<T, AllocatorDeleter<Alloc>> MakeUniqueWithAllocator(
    Args&&... args) {
  return std::unique_ptr<T, AllocatorDeleter<Alloc>>(AssertAllocatedBy<Alloc>(
      new (AllocateWith<Alloc>{}) T(std::forward<Args>(args)...)));
}

}  // namespace base
}  // namespace principia
