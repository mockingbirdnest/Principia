#pragma once

#include <memory>
#include <type_traits>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
class not_null {
 public:
  explicit not_null(Pointer const& pointer);
  explicit not_null(Pointer&& pointer);

  not_null() = delete;
  not_null(not_null const&) = default;
  // We really want a move constructor. Moving the |pointer_| might make the
  // argument an invalid |not_null<Pointer>| by making its |pointer_| null.  The
  // standard says we must leave the argument in a leave the argument in some
  // valid but otherwise indeterminate state, so we infrige on that.  On the
  // other hand, if you rely on the value after a move, you're asking for
  // trouble.
  // NOTE(egg): We would use |= default|, but VS2013 does not implement that for
  // the move constructor;
  not_null(not_null&&);
  ~not_null() = default;

  // Returns |pointer_|, by const reference to avoid issues with |unique_ptr|.
  operator Pointer const&() const;

  not_null& operator=(not_null const&) = default;
  // Implemented as a swap, so the argument remains valid.
  not_null& operator=(not_null&& other);

  // The following operators can be inlined and will be turned into constexprs
  // eventually, so having them explicitly should improve optimization.
  // They return |false|.
  bool operator==(nullptr_t other) const;
  operator bool() const;

 private:
  Pointer pointer_;
};

template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer);

// Factories taking advantage of template argument deduction.

template<typename Pointer>
not_null<Pointer> check_not_null(Pointer&& pointer);
template<typename Pointer>
not_null<Pointer> check_not_null(Pointer const& pointer);
// While the above would cover this case using the implicit conversion, this
// results in a redundant |CHECK|.
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const& pointer);
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer>&& pointer);

}  // namespace base
}  // namespace principia

#include "base/not_null_body.hpp"
