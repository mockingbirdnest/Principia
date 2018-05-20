
#pragma once

// This code comes from:
// https://code.google.com/p/or-tools/source/browse/trunk/src/base/fingerprint2011.h
// and was adapted to Visual Studio and to the needs of this project.

// Copyright 2010-2014 Google
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cstdint>

#include "base/array.hpp"

namespace principia {
namespace base {

inline std::uint64_t FingerprintCat2011(std::uint64_t const fp1,
                                        std::uint64_t const fp2) {
  // Two big prime numbers.
  std::uint64_t const mul1 = 0xC6A4A7935BD1E995u;
  std::uint64_t const mul2 = 0x228876A7198B743u;
  std::uint64_t const a = fp1 * mul1 + fp2 * mul2;
  // Note: The following line also makes sure we never return 0 or 1, because we
  // will only add something to 'a' if there are any MSBs (the remaining bits
  // after the shift) being 0, in which case wrapping around would not happen.
  return a + (~a >> 47);
}

// This should be better (collision-wise) than the default hash<std::string>,
// without being much slower. It never returns 0 or 1.
// TODO(egg): benchmark and remove the gratuitous strict aliasing violation.
inline std::uint64_t Fingerprint2011(char const* bytes, std::size_t const len) {
  // Some big prime number.
  std::uint64_t fp = 0xA5B85C5E198ED849u;
  char const* end = bytes + len;
  while (bytes + 8 <= end) {
    fp = FingerprintCat2011(
             fp, *(reinterpret_cast<const std::uint64_t*>(bytes)));
    bytes += 8;
  }
  // Note: we don't care about "consistency" (little or big endian) between
  // the bulk and the suffix of the message.
  std::uint64_t last_bytes = 0;
  while (bytes < end) {
    last_bytes += *bytes;
    last_bytes <<= 8;
    bytes++;
  }
  return FingerprintCat2011(fp, last_bytes);
}

inline std::uint64_t Fingerprint2011(Array<std::uint8_t const> const bytes) {
  return Fingerprint2011(reinterpret_cast<char const*>(bytes.data), bytes.size);
}

}  // namespace base
}  // namespace principia
