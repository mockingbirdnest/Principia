#pragma once

#include "numerics/ulp_distance.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace principia {
namespace numerics {

inline std::int64_t ULPDistance(double const x, double const y) {
  if (x == y) {
    return 0;
  }
  double const x_sign = std::copysign(1, x);
  double const y_sign = std::copysign(1, y);
  if (x_sign != y_sign) {
    double const positive = x_sign == 1 ? x : y;
    double const negative = x_sign == 1 ? y : x;
    std::int64_t const positive_distance = ULPDistance(positive, +0.0);
    std::int64_t const negative_distance = ULPDistance(negative, -0.0);
    if (positive_distance >
            std::numeric_limits<std::int64_t>::max() - negative_distance) {
      return std::numeric_limits<std::int64_t>::max();
    } else {
      return positive_distance + negative_distance;
    }
  }
  std::int64_t x_bits;
  std::int64_t y_bits;
  static_assert(sizeof(x_bits) == sizeof(x),
                "Conversion between types of different sizes");
  std::memcpy(&x_bits, &x, sizeof(x));
  static_assert(sizeof(y_bits) == sizeof(y),
                "Conversion between types of different sizes");
  std::memcpy(&y_bits, &y, sizeof(y));
  return std::abs(x_bits - y_bits);
}

}  // namespace numerics
}  // namespace principia
