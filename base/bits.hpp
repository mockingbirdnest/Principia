#pragma once

namespace principia {
namespace base {

// Floor log2 of n, or 0 for n = 0.  8 -> 3, 7 -> 2.
constexpr int FloorLog2(int n);

// Greatest power of 2 less than or equal to n.  8 -> 8, 7 -> 4.
constexpr int PowerOf2Le(int n);

}  // namespace base
}  // namespace principia

#include "base/bits_body.hpp"