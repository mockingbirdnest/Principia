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
  CHECK_LT(0, start);

  std::int64_t exponent;
  frexp(static_cast<cpp_bin_float_50>(start), &exponent);
  cpp_rational const increment =
      exp2(exponent - std::numeric_limits<double>::digits);
  cpp_rational x = start;
  // TODO(phl)Negative increment.
  for (;;) {
    auto const y = functions[0](x);
    //LOG(ERROR)<<std::setprecision(20) << y;
    std::int64_t y_exponent;
    auto const y_mantissa = frexp(y, &y_exponent);
    auto const y_scaled =
        ldexp(y_mantissa, std::numeric_limits<double>::digits);
    //LOG(ERROR)<<std::setprecision(20) <<y_scaled;
    auto const y_post_mantissa = y_scaled - floor(y_scaled);
    //LOG(ERROR)<<std::setprecision(20) << y_post_mantissa;
    auto const y_zeroes_maybe = ldexp(y_post_mantissa, zeroes);
      //LOG(ERROR)<<std::setprecision(20) << y_zeroes_maybe;
    if (trunc(y_zeroes_maybe) == 0) {
      LOG(ERROR)<<std::setprecision(20) << x << " " << y;
      return y;
    }
    x += increment;
  }
}

}  // namespace internal
}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
