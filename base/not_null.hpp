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
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, Pointer>::value>::type>
  not_null(not_null<OtherPointer> const& other);

  // We really want a move constructor. Moving the |pointer_| might make the
  // argument an invalid |not_null<Pointer>| by making its |pointer_| null.  The
  // standard says we must leave the argument in some  valid but otherwise
  // indeterminate state, so we infrige on that.  On the other hand, if you rely
  // on the value after a move, you're asking for trouble.
  // NOTE(egg): We would use |= default|, but VS2013 does not implement that for
  // the move constructor;
  not_null(not_null&&);  // NOLINT(build/c++11)
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, Pointer>::value>::type>
  not_null(not_null<OtherPointer>&& other);
  ~not_null() = default;

  not_null& operator=(not_null const&) = default;
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, Pointer>::value>::type>
  not_null& operator=(not_null<OtherPointer> const& other);
  // Implemented as a swap, so the argument remains valid.
  not_null& operator=(not_null&& other);  // NOLINT(build/c++11)
  // This operator may invalidate its argument.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, Pointer>::value>::type>
  not_null& operator=(not_null<OtherPointer>&& other);

  // Returns |pointer_|, by const reference to avoid a copy if |Pointer| is
  // |unique_ptr|.
  operator Pointer const&() const;
  // Returns |*pointer|.  Dereferencing does not work for |unique_ptr| through
  // the above implicit conversion.
  decltype(*Pointer{}) operator*() const;

  // The following operators are redundant for valid |not_null<Pointer>|s with
  // the implicit conversion to |Pointer|, but they should allow some
  // compile-time optimization.

  // Returns |false|.
  bool operator==(nullptr_t const other) const;
  // Returns |true|.
  bool operator!=(nullptr_t const other) const;
  // Returns |true|.
  operator bool() const;

 private:
  Pointer pointer_;

  template<typename OtherPointer>
  friend class not_null;
};

// Use |not_null<Pointer> const| instead.
template<typename Pointer>
class not_null<Pointer const>;

// Use |not_null<Pointer>| instead.
template<typename Pointer>
class not_null<not_null<Pointer>>;

// Use |not_null<Pointer>| instead.
template<typename Pointer>
class not_null<not_null<Pointer>&>;

// Type trait.
template<typename Pointer>
struct is_not_null : std::false_type {};
template<typename Pointer>
struct is_not_null<not_null<Pointer>> : std::true_type {};

// Factories taking advantage of template argument deduction.  They call the
// corresponding constructors for |not_null<Pointer>|.

// Returns a |not_null<Pointer>| to |*pointer|.  |CHECK|s that |pointer| is not
// null.
template<typename Pointer,
         typename = typename std::enable_if<!is_not_null<
             typename std::remove_reference<Pointer>::type>::value>::type>
not_null<typename std::remove_reference<Pointer>::type>
check_not_null(Pointer const& pointer);
// Returns a |not_null<Pointer>| to |*pointer|.  |pointer| may be invalid after
// the call, as described above.  |CHECK|s that |pointer| is not null.
template<typename Pointer,
         typename = typename std::enable_if<!is_not_null<
             typename std::remove_reference<Pointer>::type>::value>::type>
not_null<typename std::remove_reference<Pointer>::type>
check_not_null(Pointer&& pointer);  // NOLINT(build/c++11)

// While the above factories would cover this case using the implicit
// conversion, this results in a redundant |CHECK|.  These functions return
// their argument.

// Returns a copy of |pointer|.
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const& pointer);
// Equivalent to |std::move(pointer)|.  |pointer| may be invalid after the call,
// as described above.
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
