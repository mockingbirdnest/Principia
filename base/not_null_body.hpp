#pragma once

#include "base/not_null.hpp"

#include <algorithm>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
not_null<Pointer>::not_null(Pointer const& pointer) : pointer_(pointer) {}

template<typename Pointer>
not_null<Pointer>::not_null(Pointer&& pointer)  // NOLINT(build/c++11)
    : pointer_(std::move(pointer)) {}

template<typename Pointer>
template<typename OtherPointer,
         typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : pointer_(other.pointer_) {}

template<typename Pointer>
not_null<Pointer>::not_null(not_null&& other)  // NOLINT(build/c++11)
    : pointer_(std::move(other.pointer_)) {}

template<typename Pointer>
template<typename OtherPointer,
         typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)  // NOLINT
    : pointer_(std::move(other.pointer_)) {}

template<typename Pointer>
template<typename OtherPointer,
         typename>
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
template<typename OtherPointer,
         typename>
not_null<Pointer>& not_null<Pointer>::operator=(
    not_null<OtherPointer>&& other) {  // NOLINT(build/c++11)
  pointer_ = std::move(other.pointer_);
  return *this;
}

template<typename Pointer>
not_null<Pointer>::operator Pointer const&() const {
  return pointer_;
}

template<typename Pointer>
decltype(*Pointer{}) not_null<Pointer>::operator*() const {
  return *pointer_;
}

template<typename Pointer>
decltype(&*Pointer{}) const not_null<Pointer>::operator->() const {
  return &*pointer_;
}

template<typename Pointer>
template<typename P, typename>
not_null<decltype(P{}.get())> const not_null<Pointer>::get() const {
  // NOTE(egg): no |CHECK| is performed.
  return not_null<decltype(P{}.get())>(pointer_.get());
}

template<typename Pointer>
bool not_null<Pointer>::operator==(nullptr_t const other) const {
  return false;
}

template<typename Pointer>
bool not_null<Pointer>::operator!=(nullptr_t const other) const {
  return true;
}

template<typename Pointer>
not_null<Pointer>::operator bool() const {
  return true;
}

template<typename Pointer>
_checked_not_null<Pointer> check_not_null(Pointer const& pointer) {
  CHECK(pointer != nullptr);
  return not_null<typename std::remove_reference<Pointer>::type>(pointer);
}

template<typename Pointer>
_checked_not_null<Pointer>
check_not_null(Pointer&& pointer) {  // NOLINT(build/c++11)
  CHECK(pointer != nullptr);
  return not_null<typename std::remove_reference<Pointer>::type>(
      std::move(pointer));
}

template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const& pointer) {
  return pointer;
}

template<typename Pointer>
not_null<Pointer> check_not_null(
    not_null<Pointer>&& pointer) {  // NOLINT(build/c++11)
  return std::move(pointer);
}

template<typename T, typename... Args>
not_null<std::unique_ptr<T>> make_not_null_unique(Args&&... args) {  // NOLINT
  return not_null<std::unique_ptr<T>>(
      std::make_unique<T>(std::forward<Args...>(args...)));  // NOLINT
}

template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer) {
  return stream << (Pointer)(pointer);
}

}  // namespace base
}  // namespace principia
