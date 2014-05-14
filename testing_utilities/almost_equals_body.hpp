#pragma once

#include <float.h>
#include <math.h>
#include <stdint.h>

#include <string>

#include "gmock/gmock.h"

#include "quantities/dimensionless.hpp"

namespace principia {
namespace test_utilities {

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

template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected,
    int64_t const max_ulps) {
  return testing::MakePolymorphicMatcher(
      AlmostEqualsMatcher<T>(expected, max_ulps));
}

template<typename T>
AlmostEqualsMatcher<T>::AlmostEqualsMatcher(T const& expected,
                                            int64_t const max_ulps)
    : expected_(expected),
      max_ulps_(max_ulps) {}

template<typename T>
template<typename Dimensions>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return true;
  }
  int64_t const distance = ULPDistance(DoubleValue(actual),
                                       DoubleValue(expected_));
  bool const matches = distance <= max_ulps_;
  if (!matches) {
    *listener << "the numbers are separated by " << distance << " ULPs";
  }
  return matches;
}

template<typename T>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    quantities::Dimensionless const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return true;
  }
  int64_t const distance = ULPDistance(actual.value(),
                                       expected_.value());
  bool const matches = distance <= max_ulps_;
  if (!matches) {
    *listener << "the numbers are separated by " << distance << " ULPs";
  }
  return matches;;
}

template<typename T>
template<typename Scalar>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::R3Element<Scalar> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return true;
  }
  bool const x_matches =
    testing::internal::Double(DoubleValue(actual.x)).AlmostEquals(
        testing::internal::Double(DoubleValue(expected_.x)));
  bool const y_matches =
    testing::internal::Double(DoubleValue(actual.y)).AlmostEquals(
        testing::internal::Double(DoubleValue(expected_.y)));
  bool const z_matches =
    testing::internal::Double(DoubleValue(actual.z)).AlmostEquals(
        testing::internal::Double(DoubleValue(expected_.z)));
  bool const matches = x_matches && y_matches && z_matches;
  if (!matches) {
  *listener << "the following components differ: "
      << (x_matches ? "" : "x ") << (y_matches ? "" : "y ") <<
      << (z_matches ? "" : "z ");
  }
  return matches;
}

template<typename Scalar, typename Frame, unsigned int Rank>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::Multivector<Scalar, Frame, Rank> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return true;
  }
}

template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeTo(std::ostream* os) const {
  *os << "is within 4 ULPs of " << expected_;
}

template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeNegationTo(std::ostream* os) const {
  *os << "is not within 4 ULPs of " << expected_;
}

}  // namespace test_utilities
}  // namespace principia
