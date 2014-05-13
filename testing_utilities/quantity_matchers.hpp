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

}  // namespace test_utilities
}  // namespace principia
