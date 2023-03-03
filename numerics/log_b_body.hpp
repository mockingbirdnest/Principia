#pragma once

#include "numerics/log_b.hpp"

#include <limits>

namespace principia {
namespace numerics {
namespace _log_b {
namespace internal {

template<typename SourceFormat, typename>
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
    --result;
  }
  while (scaled_abs_x >= std::numeric_limits<SourceFormat>::radix) {
    scaled_abs_x /= std::numeric_limits<SourceFormat>::radix;
    ++result;
  }
  return result;
}

}  // namespace internal
}  // namespace _log_b
}  // namespace numerics
}  // namespace principia
