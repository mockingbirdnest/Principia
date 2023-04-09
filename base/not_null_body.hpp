#pragma once

#include "base/not_null.hpp"

#include <algorithm>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {
namespace _not_null {
namespace internal {

template<typename Pointer>
struct is_unique : std::false_type, not_constructible {};

template<typename T>
struct is_unique<std::unique_ptr<T>> : std::true_type, not_constructible {};

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : storage_(static_cast<pointer>(other.storage_.pointer)) {}

template<typename Pointer>
template<typename OtherPointer, typename, typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : storage_(static_cast<pointer>(other.storage_.pointer)) {}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(OtherPointer other)
    : storage_(std::move(other)) {
  CHECK(storage_.pointer != nullptr);
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)
    : storage_(std::move(other.storage_.pointer)) {}

template<typename Pointer>
template<typename OtherPointer, typename, typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)
    : storage_(static_cast<pointer>(std::move(other.storage_.pointer))) {}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null<OtherPointer> const& other) {
  storage_.pointer = other.storage_.pointer;
  return *this;
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null<OtherPointer>&& other) {
  storage_.pointer = std::move(other.storage_.pointer);
  return *this;
}

template<typename Pointer>
not_null<Pointer>::operator pointer const&&() const& {
  // This |move| is deceptive: we are not actually moving anything (|*this| is
  // |const&|), we are simply casting to an rvalue reference.
  return std::move(storage_.pointer);
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::operator OtherPointer() const& {
  return storage_.pointer;
}

template<typename Pointer>
not_null<Pointer>::operator pointer&&() && {
  return std::move(storage_.pointer);
}

template<typename Pointer>
template<typename OtherPointer, typename>
not_null<Pointer>::operator OtherPointer() && {
  return std::move(storage_.pointer);
}

template<typename Pointer>
std::add_lvalue_reference_t<typename not_null<Pointer>::element_type>
not_null<Pointer>::operator*() const {
  return *storage_.pointer;
}

template<typename Pointer>
std::add_pointer_t<typename not_null<Pointer>::element_type>
not_null<Pointer>::operator->() const {
  return std::addressof(*storage_.pointer);
}

template<typename Pointer>
template<typename P, typename>
not_null<decltype(std::declval<P>().get())> not_null<Pointer>::get() const {
  // NOTE(egg): no |CHECK| is performed.
  using type = decltype(std::declval<P>().get());
  return not_null<type>(storage_.pointer.get(), not_null<type>::unchecked_tag_);
}

template<typename Pointer>
template<typename P, typename>
not_null<decltype(std::declval<P>().release())> not_null<Pointer>::release() {
  return not_null<decltype(std::declval<P>().release())>(
      storage_.pointer.release());
}

template<typename Pointer>
template<typename Q, typename P, typename>
void not_null<Pointer>::reset(not_null<Q> const ptr) {
  storage_.pointer.reset(ptr);
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
  return storage_.pointer == other;
}

template<typename Pointer>
bool not_null<Pointer>::operator==(not_null const other) const {
  return storage_.pointer == other.storage_.pointer;
}

template<typename Pointer>
bool not_null<Pointer>::operator!=(pointer const other) const {
  return storage_.pointer != other;
}

template<typename Pointer>
bool not_null<Pointer>::operator!=(not_null const other) const {
  return storage_.pointer != other.storage_.pointer;
}

template<typename Pointer>
bool not_null<Pointer>::operator<(not_null const other) const {
  return storage_.pointer < other.storage_.pointer;
}

template<typename Pointer>
bool not_null<Pointer>::operator<=(not_null const other) const {
  return storage_.pointer <= other.storage_.pointer;
}

template<typename Pointer>
bool not_null<Pointer>::operator>=(not_null const other) const {
  return storage_.pointer >= other.storage_.pointer;
}

template<typename Pointer>
bool not_null<Pointer>::operator>(not_null const other) const {
  return storage_.pointer > other.storage_.pointer;
}

template<typename Pointer>
not_null<Pointer>::not_null(pointer other, unchecked_tag const tag)
    : storage_(std::move(other)) {}

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
  Result owned_pointer(dynamic_cast<typename Result::pointer>(unowned_pointer));
  return std::move(owned_pointer);
}

}  // namespace internal
}  // namespace _not_null
}  // namespace base
}  // namespace principia
