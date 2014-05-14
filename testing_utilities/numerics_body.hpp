#pragma once

#include <cmath>
#include <cstdint>

namespace principia {
namespace testing_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar) {
  return (scalar / Scalar::SIUnit()).value();
}

inline double RelativeError(double const expected, double const actual) {
  return std::abs(expected - actual) / std::abs(expected);
}

union Qword {
  double double_value;
  int64_t long_value;
};

int64_t ULPDistance(double const x, double const y) {
  if (x == y) {
    return 0;
  }
  if (std::copysign(x, 1) != std::copysign(y, 1)) {
    double const positive = std::copysign(x, 1) == 1 ? x : y;
    double const negative = std::copysign(x, 1) == 1 ? y : x;
    return ULPDistance(positive, +0.0) + ULPDistance(negative, -0.0);
  }
  Qword x_qword;
  Qword y_qword;
  x_qword.double_value = x;
  y_qword.double_value = y;
  return std::abs(x_qword.long_value - y_qword.long_value);
}

}  // namespace testing_utilities
}  // namespace principia
