
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

#include "base/status_or.hpp"
#include "glog/logging.h"

namespace principia {
namespace base {

template<typename T>
StatusOr<T>::StatusOr()
    : status_(Status::UNKNOWN) {}

template<typename T>
StatusOr<T>::StatusOr(Status status) : status_(std::move(status)) {
  CHECK(!status_.ok()) << "Status::OK is not a valid argument";
}

template<typename T>
StatusOr<T>::StatusOr(T const& value) {
  status_ = Status::OK;
  value_ = value;
}

template<typename T>
template<typename U>
StatusOr<T>::StatusOr(StatusOr<U> const& other)
    : status_(other.status_) {
  if (other.value_) {
    value_ = *other.value_;
  }
}

template<typename T>
template<typename U>
StatusOr<T>& StatusOr<T>::operator=(StatusOr<U> const& other) {
  status_ = other.status_;
  if (other.value_) {
    value_ = *other.value_;
  }
  return *this;
}

template<typename T>
Status const& StatusOr<T>::status() const {
  return status_;
}

template<typename T>
bool StatusOr<T>::ok() const {
  return status().ok();
}

template<typename T>
T const& StatusOr<T>::ValueOrDie() const {
  if (!status_.ok()) {
    LOG(FATAL) <<status_;
  }
  return *value_;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, StatusOr<T> const& x) {
  os << x.status();
  return os;
}

}  // namespace base
}  // namespace principia
