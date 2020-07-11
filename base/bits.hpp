
#pragma once

namespace principia {
namespace base {

// Floor log2 of n, or 0 for n = 0.  8 ↦ 3, 7 ↦ 2.
constexpr int FloorLog2(int n);

// Greatest power of 2 less than or equal to n.  8 ↦ 8, 7 ↦ 4.
constexpr int PowerOf2Le(int n);

// Computes bitreversed(bitreverse(n) + 1) assuming that n is represented on the
// given number of bits.  For 4 bits:
//   0 ↦ 8 ↦ 4 ↦ C ↦ 2 ↦ A ↦ 6 ↦ E ↦ 1 ↦ 9 ↦ 5 ↦ D ↦ 3 ↦ B ↦ 7 ↦ F
constexpr int BitReversedIncrement(int n, int bits);

}  // namespace base
}  // namespace principia

#include "base/bits_body.hpp"
