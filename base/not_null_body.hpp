#pragma once

#include "base/not_null.hpp"

#include <algorithm>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {

  /*
template<typename T>
not_null<Pointer>::not_null(not_null&&) : pointer_(std::move(pointer)) {
  CHECK(pointer_ != nullptr);
}
*/
template<typename T>
not_null<T*>::not_null(T* pointer) : pointer_(CHECK_NOTNULL(pointer)) {};

template<typename T>
not_null<T*>& not_null<T*>::operator=(not_null&& other) {
  std::swap(pointer_, other.pointer_);
  return *this;
}

template<typename T>
not_null<T*>::operator T* const() const {
  return pointer_;
}

template<typename T>
bool not_null<T*>::operator==(nullptr_t other) const {
  return false;
}

template<typename T>
not_null<T*>::operator bool() const {
  return false;
}

template<typename Pointer>
std::ostream& operator<<(std::ostream& stream, not_null<Pointer> pointer) {
  return stream << (Pointer)(pointer);
}

template<typename Pointer>
not_null<Pointer> check_not_null(Pointer pointer) {
  return not_null<Pointer>(pointer);
}

template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const pointer) {
  return pointer;
}

}  // namespace base
}  // namespace principia
