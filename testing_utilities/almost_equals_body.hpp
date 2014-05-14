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
AlmostEqualsMatcher<T>::AlmostEqualsMatcher(T expected) : expected_(expected) {}

template<typename T>
template<typename Dimensions>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> actual,
    testing::MatchResultListener * listener) const {
  bool const matches =
      testing::internal::Double(DoubleValue(actual)).AlmostEquals(
          testing::internal::Double(DoubleValue(expected)));
  if (!matches) {
    *result_listener << "the relative error is " <<
        RelativeError(DoubleValue(expected), DoubleValue(arg));
  }
  return matches;
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
