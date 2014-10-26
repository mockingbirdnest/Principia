#pragma once

#include "base/not_null.hpp"

#include <algorithm>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
not_null<Pointer>::not_null(not_null&& other)
  : pointer_(std::move(other.pointer_)) {
  CHECK(pointer_ != nullptr);
}

template<typename Pointer>
not_null<Pointer>::not_null(Pointer const& pointer) : pointer_(pointer) {
  CHECK(pointer != nullptr);
};

template<typename Pointer>
not_null<Pointer>::not_null(Pointer&& pointer) : pointer_(std::move(pointer)) {
  CHECK(pointer != nullptr);
};

template<typename Pointer>
not_null<Pointer>& not_null<Pointer>::operator=(not_null&& other) {
  std::swap(pointer_, other.pointer_);
  return *this;
}

template<typename Pointer>
not_null<Pointer>::operator Pointer const&() const {
  return pointer_;
}

template<typename Pointer>
bool not_null<Pointer>::operator==(nullptr_t other) const {
  return false;
}

template<typename Pointer>
not_null<Pointer>::operator bool() const {
  return false;
}

template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer) {
  return stream << (Pointer)(pointer);
}

template<typename Pointer>
not_null<Pointer> check_not_null(Pointer&& pointer) {
  return not_null<Pointer>(std::move(pointer));
}

template<typename Pointer>
not_null<Pointer> check_not_null(Pointer const& pointer) {
  return not_null<Pointer>(pointer);
}

template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const& pointer) {
  return pointer;
}

template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer>&& pointer) {
  return std::move(pointer)
}

}  // namespace base
}  // namespace principia
