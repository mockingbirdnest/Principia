#pragma once

#include "base/not_null.hpp"

#include <algorithm>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
not_null<Pointer>::not_null(Pointer const& pointer) : pointer_(pointer) {
  CHECK(pointer != nullptr);
}

template<typename Pointer>
not_null<Pointer>::not_null(Pointer&& pointer)  // NOLINT(build/c++11)
    : pointer_(std::move(pointer)) {
  CHECK(pointer != nullptr);
}

template<typename Pointer>
template<typename OtherPointer,
         typename>
not_null<Pointer>::not_null(not_null<OtherPointer> const& other)
    : pointer_(other.pointer_) {};

template<typename Pointer>
template<typename OtherPointer,
         typename>
not_null<Pointer>::not_null(not_null<OtherPointer>&& other)
    : pointer_(std::move(other.pointer_)) {};

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
    not_null<OtherPointer>&& other) {
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
bool not_null<Pointer>::operator==(nullptr_t const other) const {
  return false;
}

template<typename Pointer>
not_null<Pointer>::operator bool() const {
  return true;
}

template<typename Pointer>
not_null<typename remove_not_null<std::remove_reference<Pointer>::type>::type>
check_not_null(Pointer const& pointer) {
  return not_null<
      typename remove_not_null<
          std::remove_reference<Pointer>::type>::type>(pointer);
}

template<typename Pointer>
not_null<typename remove_not_null<std::remove_reference<Pointer>::type>::type>
check_not_null(Pointer&& pointer) {  // NOLINT(build/c++11)
  return not_null<
      typename remove_not_null<
          std::remove_reference<Pointer>::type>::type>(std::move(pointer));
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

template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer) {
  return stream << (Pointer)(pointer);
}

}  // namespace base
}  // namespace principia
