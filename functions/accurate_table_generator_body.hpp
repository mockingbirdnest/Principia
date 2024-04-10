#pragma once

#include "functions/accurate_table_generator.hpp"

#include "glog/logging.h"

namespace principia {
namespace functions {
namespace _accurate_table_generator {
namespace internal {

template<std::int64_t zeroes>
cpp_bin_float_50 ExhaustiveSearch(
    std::vector<AccurateFunction> const& functions,
    cpp_rational const& start) {
  std::int64_t exponent;
  frexp(static_cast<cpp_bin_float_50>(start), &exponent);
  cpp_rational const increment =
      exp2(exponent - std::numeric_limits<double>::digits);
  cpp_rational x = start;
  // TODO(phl)Negative increment.
  for (;;) {
    auto const y = functions[0](x);
    auto const y_scaled =
        ldexp(y, std::numeric_limits<double>::digits - exponent);
    auto const y_post_mantissa = y_scaled - trunc(y_scaled);
    auto const y_zeroes_maybe = ldexp(y_post_mantissa, zeroes);
    if (trunc(y) == 0) {
      LOG(ERROR) << x << " " << y;
      return y;
    }
    x += increment;
  }
}

}  // namespace internal
}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
