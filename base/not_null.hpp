#pragma once

#include <memory>
#include <type_traits>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
class not_null {
 public:
  // Creates a |not_null<Pointer>| whose |pointer_| equals the given |pointer|,
  // dawg, and checks (using a glog |CHECK| macro) that |pointer| is not null.
  explicit not_null(Pointer const& pointer);
  explicit not_null(Pointer&& pointer);  // NOLINT(build/c++11)

  not_null() = delete;
  not_null(not_null const&) = default;
  // We really want a move constructor. Moving the |pointer_| might make the
  // argument an invalid |not_null<Pointer>| by making its |pointer_| null.  The
  // standard says we must leave the argument in some  valid but otherwise
  // indeterminate state, so we infrige on that.  On the other hand, if you rely
  // on the value after a move, you're asking for trouble.
  // NOTE(egg): We would use |= default|, but VS2013 does not implement that for
  // the move constructor;
  not_null(not_null&&);  // NOLINT(build/c++11)
  ~not_null() = default;

  not_null& operator=(not_null const&) = default;
  // Implemented as a swap, so the argument remains valid.
  not_null& operator=(not_null&& other);  // NOLINT(build/c++11)

  // Returns |pointer_|, by const reference to avoid issues with |unique_ptr|.
  operator Pointer const&() const;

  // The following operators can be inlined and will be turned into constexprs
  // eventually, so having them explicitly should improve optimization.

  // Returns |false|.
  bool operator==(nullptr_t other) const;
  // Returns |true|.
  operator bool() const;

 private:
  Pointer pointer_;
};

// Factories taking advantage of template argument deduction.  They call the
// corresponding constructors for |not_null<Pointer>|.
template<typename Pointer>
not_null<Pointer> check_not_null(Pointer const& pointer);
template<typename Pointer>
not_null<Pointer> check_not_null(Pointer&& pointer);  // NOLINT(build/c++11)

// While the above would cover this case using the implicit conversion, this
// results in a redundant |CHECK|.  These return the argument.
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const& pointer);
// Returns |pointer| by moving it, which may invalidate it as described above.
template<typename Pointer>
not_null<Pointer> check_not_null(
    not_null<Pointer>&& pointer);  // NOLINT(build/c++11)

// For logging.
template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer);

}  // namespace base
}  // namespace principia

#include "base/not_null_body.hpp"
