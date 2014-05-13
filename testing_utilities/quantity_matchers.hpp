#pragma once

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
      testing::internal::Double((arg / arg.SIUnit()).Value()).AlmostEquals(
          testing::internal::Double(
              (arg_type(expected) / arg_type(expected).SIUnit()).Value()));
  if(!matches) {
    *result_listener << "the relative error is " << 
        Abs((expected - arg) / expected).Value();
  }
  return matches;
}

MATCHER_P2(Approximates, expected, expected_relative_error,
           std::string(negation ? "does not approximate " : "approximates ") +
           testing::PrintToString(expected) + " to within "
           + testing::PrintToString(expected_relative_error)) {
  double const actual_relative_error = Abs((expected - arg) / expected).Value();
  if(actual_relative_error <= expected_relative_error) {
    *result_listener << "the relative error is " <<
      actual_relative_error;
  }
  return actual_relative_error <= expected_relative_error;
}

}  // namespace test_utilities
}  // namespace principia
