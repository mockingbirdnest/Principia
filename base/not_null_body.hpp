
#pragma once

#include "base/not_null.hpp"

#include <algorithm>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
struct is_unique : std::false_type, not_constructible {};

template<typename T>
struct is_unique<std::unique_ptr<T>> : std::true_type, not_constructible {};

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : pointer_(other.pointer_) {}

template<typename Pointer>
template<typename OtherPointer, typename, typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : pointer_(static_cast<pointer>(other.pointer_)) {}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(OtherPointer other) {
  CHECK(other != nullptr);
  pointer_ = std::move(other);
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)
    : pointer_(std::move(other.pointer_)) {}

template<typename Pointer>
template<typename OtherPointer, typename, typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)
    : pointer_(static_cast<pointer>(std::move(other.pointer_))) {}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null<OtherPointer> const& other) {
  pointer_ = other.pointer_;
  return *this;
}

template<typename Pointer>
not_null<Pointer>& not_null<Pointer>::operator=(not_null&& other) {
  std::swap(pointer_, other.pointer_);
  return *this;
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null<OtherPointer>&& other) {
  pointer_ = std::move(other.pointer_);
  return *this;
}

template<typename Pointer>
not_null<Pointer>::operator pointer const&&() const& {
  // This |move| is deceptive: we are not actually moving anything (|*this| is
  // |const&|), we are simply casting to an rvalue reference.
  return std::move(pointer_);
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::operator OtherPointer() const& {
  return pointer_;
}

template<typename Pointer>
not_null<Pointer>::operator pointer&&() && {
  return std::move(pointer_);
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::operator OtherPointer() && {
  return std::move(pointer_);
}

template<typename Pointer>
std::add_lvalue_reference_t<typename not_null<Pointer>::element_type>
not_null<Pointer>::operator*() const {
  return *pointer_;
}

template<typename Pointer>
std::add_pointer_t<typename not_null<Pointer>::element_type>
not_null<Pointer>::operator->() const {
  return std::addressof(*pointer_);
}

template<typename Pointer>
template<typename P, typename>
not_null<decltype(std::declval<P>().get())> not_null<Pointer>::get() const {
  // NOTE(egg): no |CHECK| is performed.
  using type = decltype(std::declval<P>().get());
  return not_null<type>(pointer_.get(), not_null<type>::unchecked_tag_);
}

template<typename Pointer>
template<typename P, typename>
not_null<decltype(std::declval<P>().release())> not_null<Pointer>::release() {
  return not_null<decltype(std::declval<P>().release())>(pointer_.release());
}

template<typename Pointer>
template<typename Q, typename P, typename>
void not_null<Pointer>::reset(not_null<Q> const ptr) {
  pointer_.reset(ptr);
}

template<typename Pointer>
bool not_null<Pointer>::operator==(std::nullptr_t const other) const {
  return false;
}

template<typename Pointer>
bool not_null<Pointer>::operator!=(std::nullptr_t const other) const {
  return true;
}

template<typename Pointer>
not_null<Pointer>::operator bool() const {
  return true;
}

template<typename Pointer>
bool not_null<Pointer>::operator==(pointer const other) const {
  return pointer_ == other;
}

template<typename Pointer>
bool not_null<Pointer>::operator==(not_null const other) const {
  return pointer_ == other.pointer_;
}

template<typename Pointer>
bool not_null<Pointer>::operator!=(pointer const other) const {
  return pointer_ != other;
}

template<typename Pointer>
template<typename Q>
bool not_null<Pointer>::operator==(Q const other) const {
  return pointer_ == other;
}

template<typename Pointer>
template<typename Q>
bool not_null<Pointer>::operator!=(Q const other) const {
  return pointer_ != other;
}

template<typename Pointer>
bool not_null<Pointer>::operator!=(not_null const other) const {
  return pointer_ != other.pointer_;
}

template<typename Pointer>
bool not_null<Pointer>::operator<(not_null const other) const {
  return pointer_ < other.pointer_;
}

template<typename Pointer>
bool not_null<Pointer>::operator<=(not_null const other) const {
  return pointer_ <= other.pointer_;
}

template<typename Pointer>
bool not_null<Pointer>::operator>=(not_null const other) const {
  return pointer_ >= other.pointer_;
}

template<typename Pointer>
bool not_null<Pointer>::operator>(not_null const other) const {
  return pointer_ > other.pointer_;
}

template<typename Pointer>
not_null<Pointer>::not_null(pointer other, unchecked_tag const tag)
    : pointer_(std::move(other)) {}

template<typename Pointer>
_checked_not_null<Pointer> check_not_null(Pointer pointer) {
  CHECK(pointer != nullptr);
  return not_null<typename std::remove_reference<Pointer>::type>(
      std::move(pointer),
      not_null<Pointer>::unchecked_tag_);
}

#if 0
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> pointer) {
  return std::move(pointer);
}
#endif

template<typename T, typename... Args>
not_null<std::shared_ptr<T>> make_not_null_shared(Args&&... args) {
  return not_null<std::shared_ptr<T>>(
      std::make_shared<T>(std::forward<Args>(args)...),
      not_null<std::shared_ptr<T>>::unchecked_tag_);
}

template<typename T, typename... Args>
not_null<std::unique_ptr<T>> make_not_null_unique(Args&&... args) {
  return not_null<std::unique_ptr<T>>(
      std::make_unique<T>(std::forward<Args>(args)...),
      not_null<std::unique_ptr<T>>::unchecked_tag_);
}

template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer) {
  return stream << &*pointer;
}

template<typename Result, typename T>
not_null<Result> dynamic_cast_not_null(not_null<T*> const pointer) {
  static_assert(std::is_pointer<Result>::value, "|Result| should be |U*|");
  return not_null<Result>(dynamic_cast<Result>(static_cast<T*>(pointer)));
}

template<typename Result, typename T>
not_null<Result> dynamic_cast_not_null(not_null<std::unique_ptr<T>>&& pointer) {
  static_assert(is_unique<Result>::value,
                "|Result| should be |std::unique_ptr<U>|");
  T* const unowned_pointer = pointer.release();
  Result owned_pointer(
      AssertAllocatedBy<std::allocator<typename Result::element_type>>(
          dynamic_cast<typename Result::element_type*>(unowned_pointer)));
  return std::move(owned_pointer);
}

template<typename Result, typename Alloc>
not_null<Result> dynamic_cast_not_null(
    not_null<AllocatedBy<Alloc>> const pointer) {
  static_assert(std::is_pointer<Result>::value, "|Result| should be |U*|");
  using T = typename std::allocator_traits<Alloc>::value_type;
  return not_null<Result>(dynamic_cast<Result>(static_cast<T*>(pointer)));
}

}  // namespace base
}  // namespace principia
