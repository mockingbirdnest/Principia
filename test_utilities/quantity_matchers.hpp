#pragma once

#include <string>

#include "gmock/gmock.h"

namespace principia {
namespace test_utilities {

MATCHER_P(AlmostEquals, expected,
          std::string(negation ? "is not" : "is") + " almost equal to "+
          testing::PrintToString(expected)) {
  bool const matches =
      testing::internal::Double((arg / arg.SIUnit()).Value()).AlmostEquals(
          testing::internal::Double((expected / expected.SIUnit()).Value()));
  if(!matches) {
    *result_listener << "the relative error is " <<
        ToString(Abs((expected - arg) / expected)).c_str();
  }
  return matches;
}

}  // namespace test_utilities
}  // namespace principia
