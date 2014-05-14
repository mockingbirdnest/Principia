#pragma once

#include <float.h>
#include <math.h>

#include <string>

#include "gmock/gmock.h"

#include "quantities/dimensionless.hpp"

namespace principia {
namespace test_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar) {
  return (scalar / Scalar::SIUnit()).Value();
}

double RelativeError(double const expected, double const actual) {
  return std::abs(expected - actual) / std::abs(expected);
}

template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected) {
  return testing::MakePolymorphicMatcher(AlmostEqualsMatcher<T>(expected));
}

template<typename T>
AlmostEqualsMatcher<T>::AlmostEqualsMatcher(T const& expected)
    : expected_(expected) {}

template<typename T>
template<typename Dimensions>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  bool const matches =
      testing::internal::Double(DoubleValue(actual)).AlmostEquals(
          testing::internal::Double(DoubleValue(expected_)));
  if (!matches) {
    *listener << "the relative error is " <<
        RelativeError(DoubleValue(expected_), DoubleValue(actual));
  }
  return matches;
}

template<typename T>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    quantities::Dimensionless const& actual,
    testing::MatchResultListener* listener) const {
  bool const matches =
      testing::internal::Double(actual.Value()).AlmostEquals(
          testing::internal::Double(
              quantities::Dimensionless(expected_).Value()));
  if (!matches) {
    *listener << "the relative error is " <<
        RelativeError(quantities::Dimensionless(expected_).Value(),
                      actual.Value());
  }
  return matches;
}

template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeTo(std::ostream* os) const {
  *os << "is within 4 ULPs of " << expected_ << " in the uniform norm";
}
template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeNegationTo(std::ostream* os) const {
  *os << "is not within 4 ULPs of " << expected_ << " in the uniform norm";
}

}  // namespace test_utilities
}  // namespace principia
