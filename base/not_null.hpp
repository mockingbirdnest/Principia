#pragma once

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
class not_null {
 public:
  explicit not_null(Pointer pointer);

  not_null() = delete;
  not_null(not_null const&) = default;
  // Moving the |pointer_| might make the argument an invalid
  // |not_null<Pointer>| by making its |pointer_| null.  The standard says we
  // must leave the argument in a leave the argument in some valid but otherwise
  // indeterminate state.
  // not_null(not_null&&) = delete;
  not_null(not_null&&);
  ~not_null() = default;

  not_null& operator=(not_null const&) = default;
  // Implemented as a swap, so the argument remains valid.
  not_null& operator=(not_null&&);

  // Returns |pointer_|.
  operator Pointer const() const;

  // The following operators can be inlined and will be turned into constexprs
  // eventually, so having them explicitly should improve optimization.
  // They return |false|.
  bool operator==(nullptr_t other) const;
  operator bool() const;

 private:
  Pointer pointer_;
};

// Factories taking advantage of template argument deduction.

template<typename Pointer>
not_null<Pointer> check_not_null(Pointer pointer);
// While the above would cover this case using the implicit conversion, this
// results in a redundant |CHECK|.
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const pointer);

}  // namespace base
}  // namespace principia

#include "base/not_null_body.hpp"
