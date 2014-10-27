#pragma once

#include <memory>
#include <type_traits>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
class not_null;

// Type traits.
template<typename Pointer>
struct is_not_null : std::false_type {};
template<typename Pointer>
struct is_not_null<not_null<Pointer>> : std::true_type {};

template<typename Pointer>
using _checked_not_null = typename std::enable_if<
    !is_not_null<typename std::remove_reference<Pointer>::type>::value,
    not_null<typename std::remove_reference<Pointer>::type>>::type;

// |not_null<Pointer>| is a wrapper for a non-null object of type |Pointer|.
// |Pointer| should be a C-style pointer or a smart pointer.  |Pointer| must not
// be a const, reference, rvalue reference, or |not_null|.
// |not_null<Pointer>| is moveable and may be left in an invalid state when
// moved, by making its |pointer_| null.  STL mandates that we leave it in a
// valid but otherwise indeterminate state, so we infrige on that.  On the other
// hand, we're just turning one flavour of undefined behaviour into another, so
// things should be fine.
template<typename Pointer>
class not_null {
 public:

  not_null() = delete;
  not_null(not_null const&) = default;
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, Pointer>::value>::type>
  not_null(not_null<OtherPointer> const& other);

  // This constructor may invalidate its argument.
  // NOTE(egg): We would use |= default|, but VS2013 does not implement that for
  // the move constructor;
  not_null(not_null&&);  // NOLINT(build/c++11)
  // This constructor may invalidate its argument.
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
  decltype(&*Pointer{}) const operator->() const;

  template<typename P = Pointer,
           typename = decltype(P{}.get())>
  not_null<decltype(P{}.get())> const get() const;

  // The following operators are redundant for valid |not_null<Pointer>|s with
  // the implicit conversion to |Pointer|, but they should allow some
  // optimization.

  // Returns |false|.
  bool operator==(nullptr_t const other) const;
  // Returns |true|.
  bool operator!=(nullptr_t const other) const;
  // Returns |true|.
  operator bool() const;

 private:
  // Creates a |not_null<Pointer>| whose |pointer_| equals the given |pointer|,
  // dawg.  The constructor does *not* perform a null check.  Callers must
  // perform it if needed before using it.
  explicit not_null(Pointer const& pointer);
  explicit not_null(Pointer&& pointer);  // NOLINT(build/c++11)

  Pointer pointer_;

  template<typename OtherPointer>
  friend class not_null;

  template<typename P>
  friend _checked_not_null<P> check_not_null(P const& pointer);
  template<typename P>
  friend _checked_not_null<P> check_not_null(P&& pointer);
  template<typename T, typename... Args>
  friend not_null<std::unique_ptr<T>> make_not_null_unique(Args&&... args);
};

// Use |not_null<Pointer> const| instead.
template<typename Pointer>
class not_null<Pointer const>;

// Use |not_null<Pointer>| instead.
template<typename Pointer>
class not_null<not_null<Pointer>>;

// Use |not_null<Pointer>| instead.
template<typename Pointer>
class not_null<Pointer&>;

// Use |not_null<Pointer>| instead.
template<typename Pointer>
class not_null<Pointer&&>;

// Factories taking advantage of template argument deduction.  They call the
// corresponding constructors for |not_null<Pointer>|.

// Returns a |not_null<Pointer>| to |*pointer|.  |CHECK|s that |pointer| is not
// null.
template<typename Pointer>
_checked_not_null<Pointer> check_not_null(Pointer const& pointer);
// Returns a |not_null<Pointer>| to |*pointer|.  |pointer| may be invalid after
// the call.  |CHECK|s that |pointer| is not null.
template<typename Pointer>
_checked_not_null<Pointer> check_not_null(Pointer&& pointer);  // NOLINT(build/c++11)

// While the above factories would cover this case using the implicit
// conversion, this results in a redundant |CHECK|.  These functions return
// their argument.

// Returns a copy of |pointer|.
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> const& pointer);
// Returns |std::move(pointer)|.  |pointer| may be invalid after the call.
template<typename Pointer>
not_null<Pointer> check_not_null(
    not_null<Pointer>&& pointer);  // NOLINT(build/c++11)

// Factory for a |not_null<std::unique_ptr<T>>|, forwards the arguments to the
// constructor of T.  |make_not_null_unique<T>(args)| is interchangeable with
// |check_not_null(make_unique<T>(args))|, but does not perform a |CHECK|, since
// the result of |make_unique| is not null.
template<typename T, typename... Args>
not_null<std::unique_ptr<T>> make_not_null_unique(Args&&... args);

// For logging.
template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer);

}  // namespace base
}  // namespace principia

#include "base/not_null_body.hpp"
