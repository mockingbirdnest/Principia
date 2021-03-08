#pragma once

#include <memory>

#include "base/not_constructible.hpp"

namespace principia {
namespace base {

// Wrapped pointer.
// Guaranteed to be allocated by Alloc.
template<typename Alloc>
class AllocatedBy {
 public:
  using allocator = Alloc;
  using traits = std::allocator_traits<Alloc>;
  using pointer = typename traits::pointer;
  using element_type = typename traits::value_type;

  // Default constructs as nullptr.
  AllocatedBy() : ptr_(nullptr) {}

  // Implicitly convertable from AllocatedBy<B> when B converts to Alloc.
  template<typename B, typename = typename std::is_convertible<B, Alloc>>
  AllocatedBy(AllocatedBy<B> other) : ptr_(other) {}

  // Implicitly convertable from nullptr.
  AllocatedBy(std::nullptr_t) : ptr_(nullptr) {}

  // Implicitly convertable to T*.
  operator pointer() const {
    return ptr_;
  }

  // Dereference operator.
  pointer operator->() const {
    return ptr_;
  }

 private:
  AllocatedBy(pointer ptr) : ptr_(ptr) {}
  pointer ptr_;

  template<typename B>
  friend AllocatedBy<B> AssertAllocatedBy(
      typename std::allocator_traits<B>::pointer);
};

// Factory function. By calling this, you assert that |ptr| was allocated by
// |Alloc|. Typically used with allocator new.
template<typename Alloc>
AllocatedBy<Alloc> AssertAllocatedBy(
    typename std::allocator_traits<Alloc>::pointer ptr) {
  return AllocatedBy<Alloc>(ptr);
}

// Trait struct for AllocatedBy.
// Used for template metaprogramming on AllocatedBy.
template<typename T>
struct AllocatedByTraits {};

template<typename Alloc>
struct AllocatedByTraits<AllocatedBy<Alloc>> {
  using traits = std::allocator_traits<Alloc>;
};

}  // namespace base
}  // namespace principia
