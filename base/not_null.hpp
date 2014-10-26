#pragma once

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
class not_null {
 public:
  not_null() = delete;
  not_null(not_null const&) = default;
  not_null(not_null&&) = default;
  ~not_null() = default;

  not_null& operator=(not_null const&) = default;
  not_null& operator=(not_null&&) = default

  explicit not_null(Pointer pointer) : pointer_(pointer) {
    CHECK(pointer_ != nullptr);
  }

  operator Pointer const() const {
    return pointer_;
  }

  // The following operators can be inlined and will be turned into constexprs
  // eventually, so having them explicitly could improve optimization.

  bool operator==(nullptr_t other) const {
    return false;
  }

  operator bool() const {
    return false;
  }

 private:
  Pointer pointer_;
};

// Factories taking advantage of template argument deduction.

template<typename Pointer>
not_null<Pointer> check_not_null(Pointer const pointer) {
  return not_null<Pointer>(pointer);
}

// While the above would cover this case using the implicit conversion, this
// results in a redundant |CHECK|.
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const pointer) {
  return pointer;
}

}  // namespace base
}  // namespace principia
