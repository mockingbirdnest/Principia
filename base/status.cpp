
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

#include "base/status.hpp"

#include <cstdint>
#include <cstdio>
#include <ostream>
#include <string>
#include <utility>

#include "base/macros.hpp"

namespace principia {
namespace base {

inline std::string ErrorToString(Error const error) {
  switch (error) {
    case Error::OK:
      return "OK";
    case Error::CANCELLED:
      return "CANCELLED";
    case Error::UNKNOWN:
      return "UNKNOWN";
    case Error::INVALID_ARGUMENT:
      return "INVALID_ARGUMENT";
    case Error::DEADLINE_EXCEEDED:
      return "DEADLINE_EXCEEDED";
    case Error::NOT_FOUND:
      return "NOT_FOUND";
    case Error::ALREADY_EXISTS:
      return "ALREADY_EXISTS";
    case Error::PERMISSION_DENIED:
      return "PERMISSION_DENIED";
    case Error::UNAUTHENTICATED:
      return "UNAUTHENTICATED";
    case Error::RESOURCE_EXHAUSTED:
      return "RESOURCE_EXHAUSTED";
    case Error::FAILED_PRECONDITION:
      return "FAILED_PRECONDITION";
    case Error::ABORTED:
      return "ABORTED";
    case Error::OUT_OF_RANGE:
      return "OUT_OF_RANGE";
    case Error::UNIMPLEMENTED:
      return "UNIMPLEMENTED";
    case Error::INTERNAL:
      return "INTERNAL";
    case Error::UNAVAILABLE:
      return "UNAVAILABLE";
    case Error::DATA_LOSS:
      return "DATA_LOSS";
  }
  noreturn();
}

Status::Status(Error const error, std::string const& message)
    : error_(error),
      message_(error == Error::OK ? "" : message) {}

bool Status::ok() const {
  return error_ == Error::OK;
}

Error Status::error() const {
  return error_;
}

std::string const& Status::message() const {
  return message_;
}

bool Status::operator==(Status const& s) const {
  return error_ == s.error_ && message_ == s.message_;
}

bool Status::operator!=(Status const& s) const {
  return !operator==(s);
}

void Status::Update(Status const& s) {
  if (error_ == Error::OK && s.error_ != Error::OK) {
    error_ = s.error_;
    message_ = s.message_;
  }
}

std::string Status::ToString() const {
  if (error_ == Error::OK) {
    return "OK";
  } else if (message_.empty()) {
    return ErrorToString(error_);
  } else {
    return ErrorToString(error_) + ": " + message_;
  }
}

std::ostream& operator<<(std::ostream& os, Status const& s) {
  os << s.ToString();
  return os;
}

const Status Status::OK;
const Status Status::CANCELLED = Status(Error::CANCELLED, "");
const Status Status::UNKNOWN = Status(Error::UNKNOWN, "");

}  // namespace base
}  // namespace principia
