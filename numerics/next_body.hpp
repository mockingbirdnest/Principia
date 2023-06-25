#pragma once

#include "numerics/next.hpp"

#include <limits>

#include "numerics/log_b.hpp"
#include "numerics/scale_b.hpp"

namespace principia {
namespace numerics {
namespace _next {
namespace internal {

using namespace principia::numerics::_log_b;
using namespace principia::numerics::_scale_b;

template<typename SourceFormat, typename>
constexpr SourceFormat NextUp(SourceFormat const x) {
  // Handle the edge cases in the order in which they are listed in
  // IEEE 754:2008.
  if (x == -std::numeric_limits<SourceFormat>::denorm_min()) {
    return -0.0;
  }
  if (x == 0) {
    // Note that this is independent of the sign of 0, as specified by
    // IEEE 754:2008, 5.3.1.  This differs from C++ std::nextafter, which maps
    // -0 to +0.
    return std::numeric_limits<SourceFormat>::denorm_min();
  }
  if (x == std::numeric_limits<SourceFormat>::infinity()) {
    return std::numeric_limits<SourceFormat>::infinity();
  }
  if (x == -std::numeric_limits<SourceFormat>::infinity()) {
    return std::numeric_limits<SourceFormat>::lowest();
  }
  if (x != x) {
    return x;
  }
  if (x >= -std::numeric_limits<SourceFormat>::min() &&
      x < std::numeric_limits<SourceFormat>::min()) {
    // If exther x or nextUp(x) is a subnormal number, we increment by the least
    // positive magnitude.
    return x + std::numeric_limits<SourceFormat>::denorm_min();
  }
  // For a normal number, we increment by one unit in the last place of x,
  // except in the case where x is negative and nextUp(x) is in the next binade
  // (whose ULP is smaller).
  // The C++ epsilon is the ULP of one (2u).
  SourceFormat const ulp = ScaleB(
      static_cast<SourceFormat>(std::numeric_limits<SourceFormat>::epsilon()),
      static_cast<int>(LogB(x)));
  SourceFormat const signed_significand = ScaleB(x, -static_cast<int>(LogB(x)));
  if (signed_significand == -1) {
    return x + ulp / std::numeric_limits<SourceFormat>::radix;
  }
  return x + ulp;
}

template<typename SourceFormat, typename>
constexpr SourceFormat NextDown(SourceFormat const x) {
  return -NextUp(-x);
}

}  // namespace internal
}  // namespace _next
}  // namespace numerics
}  // namespace principia
