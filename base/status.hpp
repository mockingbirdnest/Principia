
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

#include <cstdint>
#include <iosfwd>
#include <string>

#include "glog/logging.h"

// TODO(phl): Many of the functions in this file should be made constexpr.

namespace principia {
namespace base {

// See https://cloud.google.com/vision/reference/rest/v1/Code for recommended
// usage of these codes.
enum class Error : std::uint64_t {
  OK                  = 0b0,
  CANCELLED           = 0b1,
  UNKNOWN             = 0b10,
  INVALID_ARGUMENT    = 0b100,
  DEADLINE_EXCEEDED   = 0b1000,
  NOT_FOUND           = 0b10000,
  ALREADY_EXISTS      = 0b100000,
  PERMISSION_DENIED   = 0b1000000,
  UNAUTHENTICATED     = 0b10000000,
  RESOURCE_EXHAUSTED  = 0b100000000,
  FAILED_PRECONDITION = 0b1000000000,
  ABORTED             = 0b10000000000,
  OUT_OF_RANGE        = 0b100000000000,
  UNIMPLEMENTED       = 0b1000000000000,
  INTERNAL            = 0b10000000000000,
  UNAVAILABLE         = 0b100000000000000,
  DATA_LOSS           = 0b1000000000000000,
};

Error operator|(Error left, Error right);
Error operator&(Error left, Error right);
Error& operator|=(Error& left, Error right);
Error& operator&=(Error& left, Error right);

std::string ErrorToString(Error error);

class Status final {
 public:
  // Creates a "successful" status.
  Status() = default;

  Status(Error error, std::string const& message);

  // Some pre-defined Status objects.
  static const Status OK;
  static const Status CANCELLED;
  static const Status UNKNOWN;

  // Accessors.
  bool ok() const;
  Error error() const;
  std::string const& message() const;

  bool operator==(Status const& s) const;
  bool operator!=(Status const& s) const;

  void Update(Status const& s);

  // Returns a combination of the error code name and message.
  std::string ToString() const;

 private:
  Error error_ = Error::OK;
  std::string message_;
};

// Prints a human-readable representation of |s| to |os|.
std::ostream& operator<<(std::ostream& os, Status const& s);

#define CHECK_OK(value) CHECK_EQ((value), ::principia::base::Status::OK)
#define EXPECT_OK(value) \
  EXPECT_THAT((value), ::principia::testing_utilities::IsOk());

#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    ::principia::base::Status const _status = (expr);                        \
    if (!_status.ok())                                                       \
      return _status;                                                        \
  } while (false)

}  // namespace base
}  // namespace principia
