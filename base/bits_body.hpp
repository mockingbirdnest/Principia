#pragma once

#include "base/bits.hpp"

#include <cstdint>

#include "base/macros.hpp"
#include "glog/logging.h"

namespace principia {
namespace base {

constexpr int FloorLog2(int const n) {
  return n == 0 ? 0 : n == 1 ? 0 : FloorLog2(n >> 1) + 1;
}

constexpr int PowerOf2Le(int const n) {
  return n == 0 ? 0 : n == 1 ? 1 : PowerOf2Le(n >> 1) << 1;
}

constexpr int BitReversedIncrement(int const n, int const bits) {
  CONSTEXPR_DCHECK(n >= 0 && n < 1 << bits);
  CONSTEXPR_DCHECK(bits > 0 && bits < 32);
  // [War03], chapter 7.1 page 105.
  std::uint32_t mask = 0x8000'0000;
  std::uint32_t x = n << (32 - bits);
  x ^= mask;
  if ((int)x >= 0) {
    do {
      mask >>= 1;
      x ^= mask;
    } while (x < mask);
  }
  return x >> (32 - bits);
}

}  // namespace base
}  // namespace principia
