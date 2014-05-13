#pragma once

#include <string>

#include "gmock/gmock.h"

namespace principia {
namespace test_utilities {

MATCHER_P2(EqualsUpToUlps, expected, ulps,
           std::string(negation ? "is not" : "is") + " within" +
           testing::PrintToString(ulps) + " ULPs of " +
           testing::PrintToString(expected)) {
  int const distance =
      testing::internal::Double::DistanceBetweenSignAndMagnitudeNumbers(
          testing::internal::Double((arg / arg.SIUnit()).Value()).bits_,
          testing::internal::Double(
              (expected / expected.SIUnit()).Value()).bits_);
  if(distance > ulps) {
    result_listener << "the numbers are separated by " << distance << " ULPs";
  }
  return distance <= ulps;
}

}  // namespace test_utilities
}  // namespace principia
