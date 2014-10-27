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

namespace {

template<typename Pointer>
using _checked_not_null = typename std::enable_if<
    !is_not_null<typename std::remove_reference<Pointer>::type>::value,
    not_null<typename std::remove_reference<Pointer>::type>>::type;

}  // namespace

// |not_null<Pointer>| is a wrapper for a non-null object of type |Pointer|.
// |Pointer| should be a C-style pointer or a smart pointer.  |Pointer| must not
// be a const, reference, rvalue reference, or |not_null|.  |not_null<Pointer>|
// is movable and may be left in an invalid state when moved, i.e., its
// |pointer_| may become null.
template<typename Pointer>
class not_null {
 public:
  not_null() = delete;
  not_null(not_null const&) = default;
  // Copy contructor for implicitly convertible pointers.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, Pointer>::value>::type>
  not_null(not_null<OtherPointer> const& other);

  // This constructor may invalidate its argument.
  // NOTE(egg): We would use |= default|, but VS2013 does not implement that for
  // the move constructor;
  not_null(not_null&&);  // NOLINT(build/c++11)
  // Move contructor for implicitly convertible pointers. This constructor may
  // invalidate its argument.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, Pointer>::value>::type>
  not_null(not_null<OtherPointer>&& other);  // NOLINT(build/c++11)
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
  not_null& operator=(not_null<OtherPointer>&& other);  // NOLINT(build/c++11)

  // Returns |pointer_|, by const reference to avoid a copy if |Pointer| is
  // |unique_ptr|.
  operator Pointer const&() const;
  // Returns |*pointer|.
  decltype(*Pointer{}) operator*() const;
  decltype(&*Pointer{}) const operator->() const;

  // When |Pointer| has a |get()| member function, this returns
  // |pointer_.get()|.
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
  friend _checked_not_null<P> check_not_null(P&& pointer);  // NOLINT
  template<typename T, typename... Args>
  friend not_null<std::unique_ptr<T>> make_not_null_unique(
      Args&&... args);  // NOLINT(build/c++11)
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
class not_null<Pointer&&>;  // NOLINT(build/c++11)

// Factories taking advantage of template argument deduction.  They call the
// corresponding constructors for |not_null<Pointer>|.

// Returns a |not_null<Pointer>| to |*pointer|.  |CHECK|s that |pointer| is not
// null.
template<typename Pointer>
_checked_not_null<Pointer> check_not_null(Pointer const& pointer);
// Returns a |not_null<Pointer>| to |*pointer|.  |pointer| may be invalid after
// the call.  |CHECK|s that |pointer| is not null.
template<typename Pointer>
_checked_not_null<Pointer> check_not_null(Pointer&& pointer);  // NOLINT

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
not_null<std::unique_ptr<T>> make_not_null_unique(
    Args&&... args);  // NOLINT(build/c++11)

// For logging.
template<typename Pointer>
std::ostream& operator<<(std::ostream& stream,
                         not_null<Pointer> const& pointer);

}  // namespace base
}  // namespace principia

#include "base/not_null_body.hpp"
