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
    for (LogBFormat k = 0; k > N; --k) {
      if (result < std::numeric_limits<SourceFormat>::radix *
                       std::numeric_limits<SourceFormat>::min()) {
        SourceFormat radix_power = 1;
        for (; k > N; --k) {
          radix_power *= std::numeric_limits<SourceFormat>::radix;
        }
        return result / radix_power;
      }
      result /= std::numeric_limits<SourceFormat>::radix;
    }
  } else {
    for (LogBFormat k = 0; k < N; ++k) {
      result *= std::numeric_limits<SourceFormat>::radix;
    }
  }
  return result;
}

}  // namespace numerics
}  // namespace principia
