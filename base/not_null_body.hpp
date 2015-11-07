#pragma once

#include "base/not_null.hpp"

#include <algorithm>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : pointer_(other.pointer_) {}

template<typename Pointer>
template<typename OtherPointer, typename, typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : pointer_(static_cast<pointer>(other.pointer_)) {}

template<typename Pointer>
not_null<Pointer>::not_null(pointer other) {
  CHECK(other != nullptr);
  pointer_ = std::move(other);
}

template<typename Pointer>
not_null<Pointer>::not_null(not_null&& other)  // NOLINT(build/c++11)
    : pointer_(std::move(other.pointer_)) {}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)  // NOLINT
    : pointer_(std::move(other.pointer_)) {}

template<typename Pointer>
template<typename OtherPointer, typename, typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)  // NOLINT
    : pointer_(static_cast<pointer>(std::move(other.pointer_))) {}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null<OtherPointer> const& other) {
  pointer_ = other.pointer_;
  return *this;
}

template<typename Pointer>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null&& other) {  // NOLINT(build/c++11)
  std::swap(pointer_, other.pointer_);
  return *this;
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null<OtherPointer>&& other) {  // NOLINT(build/c++11)
  pointer_ = std::move(other.pointer_);
  return *this;
}

template<typename Pointer>
not_null<Pointer>::operator pointer const&() const& {
  return pointer_;
}

template<typename Pointer>
not_null<Pointer>::operator pointer&&() && {
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
typename not_null<Pointer>::unchecked_tag const
    not_null<Pointer>::unchecked_tag_ = {};

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
not_null<std::unique_ptr<T>> make_not_null_unique(Args&&... args) {  // NOLINT
  return not_null<std::unique_ptr<T>>(
      std::make_unique<T>(std::forward<Args>(args)...),  // NOLINT
      not_null<std::unique_ptr<T>>::unchecked_tag_);
}

template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer) {
  return stream << &*pointer;
}

}  // namespace base
}  // namespace principia
