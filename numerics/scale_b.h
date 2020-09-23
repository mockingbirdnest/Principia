#pragma once

#include <limits>

namespace principia {
namespace numerics {

// A constexpr implementation of the IEEE 754:2008 scaleB function.
template<typename SourceFormat,
         typename LogBFormat,
         typename = std::enable_if_t<std::is_floating_point_v<SourceFormat> &&
                                     std::is_integral_v<LogBFormat>>>
constexpr SourceFormat ScaleB(SourceFormat const x, LogBFormat const N) {
  SourceFormat result = x;
  if (N < 0) {
    for (LogBFormat k = 0; k != N; --k) {
      result /= std::numeric_limits<SourceFormat>::radix;
    }
  } else {
    for (LogBFormat k = 0; k != N; ++k) {
      result *= std::numeric_limits<SourceFormat>::radix;
    }
  }
  return result;
}

static_assert(ScaleB(3.0, 2) == 0x3p2);
static_assert(ScaleB(5.0, -2) == 0x5p-2);
static_assert(ScaleB(7.0, 0) == 0x7p0);

}  // namespace numerics
}  // namespace principia
