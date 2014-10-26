#pragma once

#include <memory>
#include <type_traits>

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Pointer>
class not_null {};

template<typename T>
class not_null<T*> {
 public:
  explicit not_null(T* pointer);

  not_null() = delete;
  not_null(not_null const&) = default;
  // Moving the |pointer_| might make the argument an invalid
  // |not_null<Pointer>| by making its |pointer_| null.  The standard says we
  // must leave the argument in a leave the argument in some valid but otherwise
  // indeterminate state, so this just performs a copy.
  //not_null(not_null&&);
  ~not_null() = default;

  // Returns |pointer_|.
  operator T* const() const;

  not_null& operator=(not_null const&) = default;
  // Implemented as a swap, so the argument remains valid.
  not_null& operator=(not_null&& other);

  // The following operators can be inlined and will be turned into constexprs
  // eventually, so having them explicitly should improve optimization.
  // They return |false|.
  bool operator==(nullptr_t other) const;
  operator bool() const;

 private:
  T* pointer_;
};

template<typename T>
class not_null<std::unique_ptr<T>> {
 public:
  explicit not_null(T* pointer);

  not_null() = delete;
  not_null(not_null const&) = delete;
  // For |unique_ptr|, a move constructor is really needed.  Moving the
  // |pointer_| might make the argument an invalid |not_null<Pointer>| by making
  // its |pointer_| null.  The standard says we must leave the argument in a
  // leave the argument in some valid but otherwise indeterminate state, so we
  // infrige on that.  On the other hand, if you rely on the value after a move,
  // you're asking for trouble.
  // NOTE(egg): We would use |= default|, but VS2013 does not implement that for
  // the move constructor.
  not_null(not_null&&);
  ~not_null() = default;

  // Returns |pointer_|.
  operator T* const() const;

  not_null& operator=(not_null const&) = delete;
  // Implemented as a swap, so the argument remains valid.
  not_null& operator=(not_null&& other);

  // The following operators can be inlined and will be turned into constexprs
  // eventually, so having them explicitly should improve optimization.
  // They return |false|.
  bool operator==(nullptr_t other) const;
  operator bool() const;

 private:
  std::unique_ptr<T> pointer_;
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
