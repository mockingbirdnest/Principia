#pragma once

#include <memory>
#include <type_traits>

#include "glog/logging.h"

// This file defines a pointer wrapper |not_null| that statically ensures
// non-nullness where possible, and performs runtime checks at the point of
// conversion otherwise.
// The point is to replace cases of undefined behaviour (dereferencing a null
// pointer) by well-defined, localized, failure.
// For instance, when dereferencing a null pointer into a reference, a segfault
// will generally not occur when the pointer is dereferenced, but where the
// reference is used instead, making it hard to track where an invariant was
// violated.
// The static typing of |not_null| also optimizes away some unneeded checks:
// a function taking a |not_null| argument will not need to check its arguments,
// the caller has to provide a |not_null| pointer instead.  If the object passed
// is already a |not_null|, no check needs to be performed.
// The syntax is as follows:
//   not_null<int*> p  // non-null pointer to an |int|.
//   not_null<std::unique_ptr<int>> p // non-null unique pointer to an |int|.
// |not_null| does not have a default constructor, since there is no non-null
// default valid pointer.  The only ways to construct a |not_null| pointer,
// other than from existing instances of |not_null|, are |check_not_null| and
// |make_unique_not_null|.
//
// The following example shows uses of |not_null|:
//   void Accumulate(not_null<int*> const accumulator,
//                   not_null<int const*> const term) {
//     *accumulator += *term;  // This will not dereference a null pointer.
//   }
//
//   void InterfaceAccumulator(int* const dubious_accumulator,
//                             int const* const term_of_dubious_c_provenance) {
//     // The call below performs two checks.  If either parameter is null, the
//     // program will fail (through a glog |CHECK|) at the callsite.
//     Accumulate(check_not_null(dubious_accumulator),
//                check_not_null(term_of_dubious_c_provenance));
//     // The call below fails to compile: we need to check the arguments.
//     Accumulate(dubious_accumulator, term_of_dubious_c_provenance);
//   }
//
//   void UseAccumulator() {
//     not_null<std::unique_ptr<int>> accumulator =
//         make_not_null_unique<int>(0);
//     not_null<int> term = // ...
//     // This compiles, no check is performed.
//     Accumulate(accumulator.get(), term);
//     // ...
//   }
//
// The following redundant checks are not performed.  This is useful, since
// in a template we may not know whether a pointer is |not_null|.
//   not_null<std::unique_ptr<int>> accumulator = // ...
//   not_null<int> term = // ...
// |check_not_null| does not perform a check.
//   Accumulate(check_not_null(accumulator.get()), check_not_null(term));
// |term == nullptr| can be expanded to false through inlining, so the branch
// will likely be optimized away.
//   if (term == nullptr) // ...
//   // Same as above.
//   if (term) // ...

namespace principia {
namespace base {

template<typename Pointer>
class not_null;

// Type traits.

// |is_instance<T, U>::value| is true if and only if |U| is an instance of the
// template |T|.  It is false otherwise.
template<template<typename...> class T, typename U>
struct is_instance_of : std::false_type {};
template<template<typename...> class T, typename U>
struct is_instance_of<T, T<U>> : std::true_type {};

// |remove_not_null<not_null<T>>::type| is |remove_not_null<T>::type|.
// The recurrence ends when |T| is not an instance of |not_null|, in which case
// |remove_not_null<T>::type| is |T|.
template<typename Pointer>
struct remove_not_null {
  using type = Pointer;
};
template<typename Pointer>
struct remove_not_null<not_null<Pointer>> {
  using type = typename remove_not_null<Pointer>::type;
};

// When |T| is not a reference, |_checked_not_null<T>| is |not_null<T>| if |T|
// is not already an instance of |not_null|.  It fails otherwise.
// |_checked_not_null| is invariant under application of reference or rvalue
// reference to its template argument.
template<typename Pointer>
using _checked_not_null = typename std::enable_if<
    !is_instance_of<not_null,
                    typename std::remove_reference<Pointer>::type>::value,
    not_null<typename std::remove_reference<Pointer>::type>>::type;

// |not_null<Pointer>| is a wrapper for a non-null object of type |Pointer|.
// |Pointer| should be a C-style pointer or a smart pointer.  |Pointer| must not
// be a const, reference, rvalue reference, or |not_null|.  |not_null<Pointer>|
// is movable and may be left in an invalid state when moved, i.e., its
// |pointer_| may become null.
// |not_null<not_null<Pointer>>| and |not_null<Pointer>| are equivalent.
// This is useful when a |template<typename T>| using a |not_null<T>| is
// instanced with an instance of |not_null|.
template<typename Pointer>
class not_null {
 public:
  // The type of the pointer being wrapped.
  // This follows the naming convention from |std::unique_ptr|.
  using pointer = typename remove_not_null<Pointer>::type;

  not_null() = delete;

  // Copy constructor from an other |not_null<Pointer>|.
  not_null(not_null const&) = default;
  // Copy contructor for implicitly convertible pointers.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, pointer>::value>::type>
  not_null(not_null<OtherPointer> const& other);
  // Copy constructor from a nullable pointer, performs a null check.
  not_null(pointer other);
  // Explicit copy constructor for static_cast'ing.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               !std::is_convertible<OtherPointer, pointer>::value>::type,
           typename = decltype(static_cast<pointer>(
                                   std::declval<OtherPointer>()))>
  explicit not_null(not_null<OtherPointer> const& other);

  // Move constructor from an other |not_null<Pointer>|.  This constructor may
  // invalidate its argument.
  // NOTE(egg): We would use |= default|, but VS2013 does not implement that for
  // the move constructor.
  not_null(not_null&&);  // NOLINT(build/c++11)
  // Move contructor for implicitly convertible pointers. This constructor may
  // invalidate its argument.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, pointer>::value>::type>
  not_null(not_null<OtherPointer>&& other);  // NOLINT(build/c++11)
  // Explicit move constructor for static_cast'ing. This constructor may
  // invalidate its argument.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               !std::is_convertible<OtherPointer, pointer>::value>::type,
           typename = decltype(static_cast<pointer>(
                                   std::declval<OtherPointer>()))>
  explicit not_null(not_null<OtherPointer>&& other);  // NOLINT(build/c++11)

  ~not_null() = default;

  // Copy assigment operators.
  not_null& operator=(not_null const&) = default;
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, pointer>::value>::type>
  not_null& operator=(not_null<OtherPointer> const& other);

  // Move assignment operators.
  // Implemented as a swap, so the argument remains valid.
  not_null& operator=(not_null&& other);  // NOLINT(build/c++11)
  // This operator may invalidate its argument.
  template<typename OtherPointer,
           typename = typename std::enable_if<
               std::is_convertible<OtherPointer, pointer>::value>::type>
  not_null& operator=(not_null<OtherPointer>&& other);  // NOLINT(build/c++11)

  // Returns |pointer_|, by const reference to avoid a copy if |pointer| is
  // |unique_ptr|.
  operator pointer const&() const;
  // Returns |*pointer_|.
  decltype(*std::declval<pointer>()) operator*() const;
  decltype(std::addressof(*std::declval<pointer>())) const operator->() const;

  // When |pointer| has a |get()| member function, this returns
  // |pointer_.get()|.
  template<typename P = pointer, typename = decltype(std::declval<P>().get())>
  not_null<decltype(std::declval<P>().get())> get() const;

  // When |pointer| has a |release()| member function, this returns
  // |pointer_.release()|.  May invalidate its argument.
  template<typename P = pointer,
           typename = decltype(std::declval<P>().release())>
  not_null<decltype(std::declval<P>().release())> release();

  // When |pointer| has a |reset()| member function, this calls
  // |pointer_.reset()|.
  template<typename Q,
           typename P = pointer,
           typename = decltype(std::declval<P>().reset())>
  void reset(not_null<Q> const ptr);

  // The following operators are redundant for valid |not_null<Pointer>|s with
  // the implicit conversion to |pointer|, but they should allow some
  // optimization.

  // Returns |false|.
  bool operator==(std::nullptr_t const other) const;
  // Returns |true|.
  bool operator!=(std::nullptr_t const other) const;
  // Returns |true|.
  operator bool() const;

  // Equality.
  bool operator==(pointer const other) const;
  bool operator==(not_null const other) const;
  bool operator!=(pointer const other) const;
  bool operator!=(not_null const other) const;

  // Ordering.
  bool operator<(not_null const other) const;
  bool operator<=(not_null const other) const;
  bool operator>=(not_null const other) const;
  bool operator>(not_null const other) const;

 private:
  struct unchecked_tag {
    inline unchecked_tag() {}
  };

  // Creates a |not_null<Pointer>| whose |pointer_| equals the given |pointer|,
  // dawg.  The constructor does *not* perform a null check.  Callers must
  // perform one if needed before using it.
  explicit not_null(pointer ptr, unchecked_tag const tag);

  pointer pointer_;

  static unchecked_tag const unchecked_tag_;

  template<typename OtherPointer>
  friend class not_null;

  template<typename P>
  friend _checked_not_null<P> check_not_null(P pointer);
  template<typename T, typename... Args>
  friend not_null<std::unique_ptr<T>> make_not_null_unique(
      Args&&... args);  // NOLINT(build/c++11)
  template<typename Pointer>
  friend std::ostream& operator<<(std::ostream& stream,
                                  not_null<Pointer> const& pointer);
};

// We want only one way of doing things, and we can't make
// |not_null<Pointer> const| and |not_null<Pointer const>| etc. equivalent
// easily.

// Use |not_null<Pointer> const| instead.
template<typename Pointer>
class not_null<Pointer const>;
// Use |not_null<Pointer>&| instead.
template<typename Pointer>
class not_null<Pointer&>;
// Use |not_null<Pointer>&&| instead.
template<typename Pointer>
class not_null<Pointer&&>;  // NOLINT(build/c++11)

// Factory taking advantage of template argument deduction.  Returns a
// |not_null<Pointer>| to |*pointer|.  |CHECK|s that |pointer| is not null.
template<typename Pointer>
_checked_not_null<Pointer> check_not_null(Pointer pointer);

#if 0
// While the above factory would cover this case using the implicit
// conversion, this results in a redundant |CHECK|.
// This function returns its argument.
template<typename Pointer>
not_null<Pointer> check_not_null(not_null<Pointer> pointer);
#endif

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
