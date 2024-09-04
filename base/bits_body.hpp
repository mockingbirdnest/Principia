#pragma once

#include "base/bits.hpp"

#include "base/macros.hpp"  // 🧙 For CONSTEXPR_DCHECK.
#include "glog/logging.h"

namespace principia {
namespace base {
namespace _bits {
namespace internal {

constexpr std::int64_t FloorLog2(std::int64_t const n) {
  return n == 0 ? 0 : n == 1 ? 0 : FloorLog2(n >> 1) + 1;
}

constexpr std::int64_t PowerOf2Le(std::int64_t const n) {
  return n == 0 ? 0 : n == 1 ? 1 : PowerOf2Le(n >> 1) << 1;
}

constexpr std::int64_t BitReversedIncrement(std::int64_t const n,
                                            std::int64_t const bits) {
  if (bits == 0) {
    CONSTEXPR_DCHECK(n == 0);
    return 0;
  }
  CONSTEXPR_DCHECK(n >= 0 && n < 1LL << bits);
  CONSTEXPR_DCHECK(bits > 0 && bits < 32);
  // [War03], chapter 7.1 page 105.
  std::uint32_t mask = 0x8000'0000;
  std::uint32_t x = n << (32 - bits);
  x ^= mask;
  if (static_cast<std::int32_t>(x) >= 0) {
    do {
      mask >>= 1;
      x ^= mask;
    } while (x < mask);
  }
  return x >> (32 - bits);
}

}  // namespace internal
}  // namespace _bits
}  // namespace base
}  // namespace principia
