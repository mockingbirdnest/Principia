#pragma once

#include "functions/accurate_table_generator.hpp"

#include "glog/logging.h"

namespace principia {
namespace functions {
namespace _accurate_table_generator {
namespace internal {

template<std::int64_t zeroes>
bool HasDesiredZeroes(cpp_bin_float_50 const& y) {
  std::int64_t y_exponent;
  auto const y_mantissa = frexp(y, &y_exponent);
  auto const y_mantissa_scaled =
      ldexp(y_mantissa, std::numeric_limits<double>::digits);
  auto const y_post_mantissa = y_mantissa_scaled - floor(y_mantissa_scaled);
  auto const y_candidate_zeroes = ldexp(y_post_mantissa, zeroes);
  return trunc(y_candidate_zeroes) == 0;
}

template<std::int64_t zeroes>
cpp_rational ExhaustiveSearch(AccurateFunction const& function,
                              cpp_rational const& start) {
  CHECK_LT(0, start);

  // We will look for candidates both above and below |start|.  Note that if
  // start is a power of 2, the increments above and below |start| are not the
  // same.
  std::int64_t exponent;
  auto const start_mantissa =
      frexp(static_cast<cpp_bin_float_50>(start), &exponent);
  cpp_rational const high_increment =
      exp2(exponent - std::numeric_limits<double>::digits);
  cpp_rational const low_increment = start_mantissa ==
      0.5 ? high_increment / 2 : high_increment;

  cpp_rational high_x = start;
  cpp_rational low_x = start - low_increment;
  for (;;) {
    if (HasDesiredZeroes<zeroes>(function(high_x))) {
      return high_x;
    }
    high_x += high_increment;
    if (HasDesiredZeroes<zeroes>(function(low_x))) {
      return low_x;
    }
    low_x -= low_increment;
  }
}

}  // namespace internal
}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
