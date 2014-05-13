#pragma once

#include <float.h>

#include <string>

#include "gmock/gmock.h"

#include "quantities/dimensionless.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace test_utilities {

MATCHER_P(AlmostEquals, expected,
          std::string(negation ? "is not" : "is") + " within 4 ULPs of " +
          testing::PrintToString(expected)) {
  bool const matches =
      testing::internal::Double(DoubleValue(arg)).AlmostEquals(
          testing::internal::Double(DoubleValue(arg_type(expected))));
  if (!matches) {
    *result_listener << "the relative error is " <<
        RelativeError(DoubleValue(arg_type(expected)), DoubleValue(arg));
  }
  return matches;
}

MATCHER_P(AlmostVanishesBefore, input_magnitude,
          std::string(negation ? "is not" : "is") + " within " +
          testing::PrintToString(input_magnitude) + " * 4 epsilon of zero.") {
  double const expected_absolute_error = input_magnitude * 4 * DBL_EPSILON;
  double const actual_absolute_error = std::abs(DoubleValue(arg));
  if (actual_absolute_error <= expected_absolute_error) {
    *result_listener << "the absolute error is " <<
      actual_absolute_error / (4 * DBL_EPSILON) << " * 4 epsilon";
  }
  return actual_absolute_error <= expected_absolute_error;
}

MATCHER_P2(Approximates, expected, expected_relative_error,
           std::string(negation ? "does not approximate " : "approximates ") +
           testing::PrintToString(expected) + " to within "
           + testing::PrintToString(expected_relative_error)) {
  double const actual_relative_error =
      RelativeError(DoubleValue(arg_type(expected)), DoubleValue(arg));
  if (actual_relative_error <= expected_relative_error) {
    *result_listener << "the relative error is " << actual_relative_error;
  }
  return actual_relative_error <= expected_relative_error;
}

template<typename Scalar>
double DoubleValue(Scalar const& scalar) {
  return (scalar / Scalar::SIUnit()).Value();
}

double RelativeError(double const expected, double const actual) {
  return std::abs(expected - actual) / std::abs(expected);
}

}  // namespace test_utilities
}  // namespace principia
