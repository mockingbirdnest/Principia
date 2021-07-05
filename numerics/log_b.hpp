#pragma once

#include <limits>

namespace principia {
namespace numerics {

// A constexpr implementation of the IEEE 754:2008 logB function.
// Uses sourceFormat as logBFormat, which makes it easy to cleanly handle NaN,
// infinity, and 0.
template<typename SourceFormat,
         typename = std::enable_if_t<std::is_floating_point_v<SourceFormat>>>
constexpr SourceFormat LogB(SourceFormat const x) {
  if (x != x) {
    return x;
  }
  SourceFormat const abs_x = x > 0 ? x : 0 - x;
  if (abs_x == std::numeric_limits<SourceFormat>::infinity()) {
    return std::numeric_limits<SourceFormat>::infinity();
  }
  if (x == 0) {
    // Signal the divideByZero exception (by dividing by -0).
    return 1 / -abs_x;
  }
  SourceFormat scaled_abs_x = abs_x;
  SourceFormat result = 0;
  while (scaled_abs_x < 1) {
    scaled_abs_x *= std::numeric_limits<SourceFormat>::radix;
    result -= 1;
  }
  while (scaled_abs_x >= std::numeric_limits<SourceFormat>::radix) {
    scaled_abs_x /= std::numeric_limits<SourceFormat>::radix;
    result += 1;
  }
  return result;
}

}  // namespace numerics
}  // namespace principia
