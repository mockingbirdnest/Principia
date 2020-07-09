#pragma once

#include "base/bits.hpp"

namespace principia {
namespace base {

constexpr int FloorLog2(int const n) {
  return n == 0 ? 0 : n == 1 ? 0 : FloorLog2(n >> 1) + 1;
}

constexpr int PowerOf2Le(int const n) {
  return n == 0 ? 0 : n == 1 ? 1 : PowerOf2Le(n >> 1) << 1;
}

}  // namespace base
}  // namespace principia
