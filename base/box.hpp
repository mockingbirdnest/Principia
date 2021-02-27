// Allocator-based unique pointer implementation.

#pragma once

#include <cstddef>
#include <memory>
#include <type_traits>

namespace principia {
namespace base {

template <typename T, typename Allocator>
class NullableBox;

// Smart pointer that uses an allocator  (as opposed to a deleter like
// unique_ptr does). The allocator has the same semantics as STL container
// allocators.
//
// Can't be null unless moved out of (same semantics as base::not_null). Do not
// use moved-from objects.
template <typename T, typename Allocator = std::allocator<T>>
class Box final {
 public:
  using pointer = T*;
  using element_type = T;
  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<Allocator>;

  // Allocates memory but does not initialize the memory. Dangerous! Used with
  // private T constructors.
  static Box<T, Allocator> UnsafeMakeUninitialized();

  // Constructs a Box from a NullableBox. CHECK-fails if null.
  static Box<T, Allocator> FromNullableChecked(
      NullableBox<T, Allocator>&& nullable);

  // Primary constructor. Forwards arguments to T's constructor.
  template <typename... Args>
  explicit Box(Args&&... args);

  // Move constructor. Nullifies |other|!
  Box(Box&& other);

  // Templated move constructor. Useful for upcasting.
  //
  // Only allowed when U* is implicitly convertible to T* and when the
  // allocators are compatible.
  template <typename U, typename E>
  Box(Box<U, E>&& other);

  // Move assignment operator. Equivalent to swap(|other|).
  Box& operator=(Box&& other);

  // Consumes this object and covnerts it into a NullableBox.
  NullableBox<T, Allocator> IntoNullable() &&;

  // Swaps with |other|.
  void swap(Box& other);

  // Returns a pointer to the managed object.
  T* get() const;

  // Dereferences pointer to the managed object.
  typename std::add_lvalue_reference<T>::type operator*() const;
  T* operator->() const;

  // Equality.
  bool operator==(T* other) const;
  bool operator==(Box const& other) const;
  bool operator!=(T* other) const;
  bool operator!=(Box const& other) const;

  // Ordering.
  bool operator<(Box const& other) const;
  bool operator<=(Box const& other) const;
  bool operator>=(Box const& other) const;
  bool operator>(Box const& other) const;

 private:
  struct UninitTag {};
  struct NullTag {};
  explicit Box(UninitTag);
  explicit Box(NullTag);

  // Invariant: not null.
  NullableBox<T, Allocator> box_;

  template <typename U, typename E>
  friend class Box;
};

template <typename T, typename Allocator = std::allocator<T>>
class NullableBox final {
 public:
  using pointer = T*;
  using element_type = T;
  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<Allocator>;

  // Constraints on Allocator type.
  static_assert(allocator_traits::is_always_equal::value,
                "Allocator must be stateless");
  static_assert(
      std::is_same<typename allocator_traits::pointer, pointer>::value,
      "Allocator pointer type must match");

  // Construct a NullableBox from a Box by throwing away non-null guarantee.
  NullableBox(Box<T>&& box);

  // Construct a NullableBox from nullptr.
  NullableBox(std::nullptr_t);

  // Move constructor. Nullifies |other|.
  NullableBox(NullableBox&& other);

  // Templated move constructor. Useful for upcasting.
  //
  // Only allowed when U* is implicitly convertible to T* and when the
  // allocators are compatible.
  template <
      typename U, typename E,
      typename = typename std::enable_if<std::is_same<
          typename allocator_traits::template rebind_alloc<U>, E>::value>::type>
  NullableBox(NullableBox<U, E>&& other);

  // Move assignment operator. Equivalent to swap(|other|).
  NullableBox& operator=(NullableBox&& other);

  // This is a move-only type.
  NullableBox(NullableBox const& other) = delete;
  NullableBox& operator=(NullableBox const& other) = delete;

  // Destructor.
  ~NullableBox();

  // Consumes this object and convert it into a Box. CHECK-fails if null.
  Box<T, Allocator> AssertNotNull() &&;

  // Swaps with |other|.
  void swap(NullableBox& other);

  // Returns a pointer to the managed object.
  T* get() const;

  // Dereferences pointer to the managed object.
  typename std::add_lvalue_reference<T>::type operator*() const;
  T* operator->() const;

  // Equality.
  bool operator==(T* other) const;
  bool operator==(NullableBox const& other) const;
  bool operator!=(T* other) const;
  bool operator!=(NullableBox const& other) const;

  // Ordering.
  bool operator<(NullableBox const& other) const;
  bool operator<=(NullableBox const& other) const;
  bool operator>=(NullableBox const& other) const;
  bool operator>(NullableBox const& other) const;

 private:
  // Allocates a T without constructing it.
  explicit NullableBox();

  // Invariant: ptr_ was allocated with Allocator or is null.
  T* ptr_;

  template <typename U, typename E>
  friend class NullableBox;

  template <typename U, typename E>
  friend class Box;
};

template <typename T, typename Allocator>
void swap(Box<T, Allocator>& lhs, Box<T, Allocator>& rhs);

template <typename T, typename Allocator>
void swap(NullableBox<T, Allocator>& lhs, NullableBox<T, Allocator>& rhs);

}  // namespace base
}  // namespace principia

#include "base/box_body.hpp"
