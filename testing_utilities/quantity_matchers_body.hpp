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

template<typename Scalar>
AlmostEqualsMatcher<Scalar>::AlmostEqualsMatcher(Scalar expected)
    : expected_(expected) {}

template<typename Scalar>
bool AlmostEqualsMatcher<Scalar>::MatchAndExplain(
    Scalar actual,
    testing::MatchResultListener * listener) const {
}

template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeTo(std::ostream* os) const;
template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeNegationTo(std::ostream* os) const;

}  // namespace test_utilities
}  // namespace principia
