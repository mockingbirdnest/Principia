#pragma once

#include "numerics/scale_b.hpp"

#include <limits>

namespace principia {
namespace numerics {
namespace _scale_b {
namespace internal {

template<typename SourceFormat, typename LogBFormat, typename>
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

}  // namespace internal
}  // namespace _scale_b
}  // namespace numerics
}  // namespace principia
