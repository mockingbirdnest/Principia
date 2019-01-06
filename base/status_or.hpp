
#pragma once

// This code comes from:
// https://github.com/google/protobuf/tree/master/src/google/protobuf/stubs
// and was adapted to Visual Studio and to the needs of this project.

// Protocol Buffers - Google's data interchange format
// Copyright 2008 Google Inc.  All rights reserved.
// https://developers.google.com/protocol-buffers/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// StatusOr<T> is the union of a Status object and a T
// object. StatusOr models the concept of an object that is either a
// usable value, or an error Status explaining why such a value is
// not present. To this end, StatusOr<T> does not allow its Status
// value to be Status::OK.
//
// The primary use-case for StatusOr<T> is as the return value of a
// function which may fail.
//
// Example client usage for a StatusOr<T>, where T is not a pointer:
//
//  StatusOr<float> result = DoBigCalculationThatCouldFail();
//  if (result.ok()) {
//    float answer = result.ValueOrDie();
//    printf("Big calculation yielded: %f", answer);
//  } else {
//    LOG(ERROR) << result.status();
//  }
//
// Example client usage for a StatusOr<T*>:
//
//  StatusOr<Foo*> result = FooFactory::MakeNewFoo(arg);
//  if (result.ok()) {
//    std::unique_ptr<Foo> foo(result.ValueOrDie());
//    foo->DoSomethingCool();
//  } else {
//    LOG(ERROR) << result.status();
//  }
//
// Example client usage for a StatusOr<std::unique_ptr<T>>:
//
//  StatusOr<std::unique_ptr<Foo>> result = FooFactory::MakeNewFoo(arg);
//  if (result.ok()) {
//    std::unique_ptr<Foo> foo = result.ConsumeValueOrDie();
//    foo->DoSomethingCool();
//  } else {
//    LOG(ERROR) << result.status();
//  }
//
// Example factory implementation returning StatusOr<T*>:
//
//  StatusOr<Foo*> FooFactory::MakeNewFoo(int arg) {
//    if (arg <= 0) {
//      return ::util::Status(::util::error::INVALID_ARGUMENT,
//                            "Arg must be positive");
//    } else {
//      return new Foo(arg);
//    }
//  }
//

#include <new>
#include <optional>
#include <string>
#include <utility>

#include "base/status.hpp"

namespace principia {
namespace base {

template<typename T>
class StatusOr final {
 public:
  // Construct a new StatusOr with Status::UNKNOWN status
  StatusOr();

  // Construct a new object with the given non-ok status. After calling
  // this constructor, calls to |ValueOrDie()| will fail.
  //
  // NOTE: Not explicit - we want to use |StatusOr<T>| as a return
  // value, so it is convenient and sensible to be able to do |return Status();|
  // when the return type is |StatusOr<T>|.
  StatusOr(Status status);  // NOLINT

  // Construct a new object with the given value. After calling this
  // constructor, calls to |ValueOrDie()| will succeed, and calls to |status()|
  // will return |Status::OK|.
  //
  // NOTE: Not explicit - we want to use StatusOr<T> as a return type
  // so it is convenient and sensible to be able to do |return T()|
  // when when the return type is |StatusOr<T>|.
  StatusOr(T const& value);  // NOLINT

  // Copy constructor.
  StatusOr(StatusOr const& other) = default;

  // Conversion copy constructor, |T| must be copy-constructible from |U|.
  template<typename U>
  StatusOr(StatusOr<U> const& other);

  // Assignment operator.
  StatusOr& operator=(StatusOr const& other) = default;

  // Conversion assignment operator, |T| must be copy-assignable from |U|.
  template<typename U>
  StatusOr& operator=(StatusOr<U> const& other);

  // Returns a reference to our status. If this contains a |T|, then
  // returns |Status::OK|.
  Status const& status() const;

  // Returns true iff the status is OK (and |ValueOrDie()| may be called).
  bool ok() const;

  // Returns a reference to our current value, or fails if the status is not OK.
  T const& ValueOrDie() const;

 private:
  Status status_;
  std::optional<T> value_;

  template<typename U>
  friend class StatusOr;
};

// Prints a human-readable representation of |x| to |os|.
template<typename T>
std::ostream& operator<<(std::ostream& os, StatusOr<T> const& x);

}  // namespace base
}  // namespace principia

#include "base/status_or_body.hpp"
