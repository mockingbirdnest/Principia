#pragma once

// This code comes from:
// https://code.google.com/p/or-tools/source/browse/trunk/src/base/fingerprint2011.h
// and was adapted to Visual Studio 2013.

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

inline std::uint64_t FingerprintCat2011(std::uint64_t fp1, std::uint64_t fp2) {
  // Two big prime numbers.
  const std::uint64_t kMul1 = 0xc6a4a7935bd1e995ULL;
  const std::uint64_t kMul2 = 0x228876a7198b743ULL;
  std::uint64_t a = fp1 * kMul1 + fp2 * kMul2;
  // Note: The following line also makes sure we never return 0 or 1, because we
  // will only add something to 'a' if there are any MSBs (the remaining bits
  // after the shift) being 0, in which case wrapping around would not happen.
  return a + (~a >> 47);
}

// This should be better (collision-wise) than the default hash<std::string>,
// without being much slower. It never returns 0 or 1.
inline std::uint64_t Fingerprint2011(const char* bytes, size_t len) {
  // Some big prime numer.
  std::uint64_t fp = 0xa5b85c5e198ed849ULL;
  const char* end = bytes + len;
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
