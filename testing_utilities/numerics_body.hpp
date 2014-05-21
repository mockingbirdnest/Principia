#pragma once

#include <cmath>
#include <cstdint>

namespace principia {
namespace testing_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar) {
  return (scalar / Scalar::SIUnit()).value();
}

template<typename T, typename Norm>
quantities::Dimensionless RelativeError(T const& expected, T const& actual,
                                        Norm const norm) {
  return norm(expected - actual) / norm(expected);
}

inline quantities::Dimensionless RelativeError(
    quantities::Dimensionless const& expected,
    quantities::Dimensionless const& actual) {
  return RelativeError(expected, actual, quantities::Abs);
}

template<typename Dimensions>
quantities::Dimensionless RelativeError(
    quantities::Quantity<Dimensions> const& expected,
    quantities::Quantity<Dimensions> const& actual) {
  return RelativeError(expected, actual, quantities::Abs<Dimensions>);
}

template<typename Scalar>
geometry::R3Element<Scalar> RelativeError(
    geometry::R3Element<Scalar> const& expected,
    geometry::R3Element<Scalar> const& actual) {
  return RelativeError(expected, actual,
                       [](R3Element<Scalar> v) { return v.Norm() });
}

union Qword {
  double double_value;
  std::int64_t long_value;
};

inline std::int64_t ULPDistance(double const x, double const y) {
  if (x == y) {
    return 0;
  }
  double const x_sign = std::copysign(1, x);
  double const y_sign = std::copysign(1, y);
  if (x_sign != y_sign) {
    double const positive = x_sign == 1 ? x : y;
    double const negative = x_sign == 1 ? y : x;
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
