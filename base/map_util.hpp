#pragma once

// This code comes from:
// https://github.com/google/protobuf/blob/master/src/google/protobuf/stubs/map_util.h
// and was adapted to Visual Studio 2013 and to the needs of this project.

// Protocol Buffers - Google's data interchange format
// Copyright 2014 Google Inc.  All rights reserved.
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

// from google3/util/gtl/map_util.h
// Author: Anton Carver

#include "glog/logging.h"

namespace principia {
namespace base {
namespace _map_util {
namespace internal {

// Returns true if and only if the given collection contains the given key.
template <class Collection>
bool Contains(Collection const& collection,
              typename Collection::key_type const key) {
  return collection.find(key) != collection.end();
}

template<class Collection>
const typename Collection::value_type::second_type& FindOrDie(
    Collection const& collection,
    typename Collection::value_type::first_type const& key) {
  typename Collection::const_iterator it = collection.find(key);
  CHECK(it != collection.end()) << "Map key not found: " << key;
  return it->second;
}

// Same as above, but returns a non-const reference.
template<class Collection>
typename Collection::value_type::second_type& FindOrDie(
    Collection& collection,
    typename Collection::value_type::first_type const& key) {
  typename Collection::iterator it = collection.find(key);
  CHECK(it != collection.end()) << "Map key not found: " << key;
  return it->second;
}

}  // namespace internal

using internal::Contains;
using internal::FindOrDie;

}  // namespace _map_util
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_map_util;
}  // namespace principia::base
