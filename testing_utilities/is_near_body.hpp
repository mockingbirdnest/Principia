
#pragma once

#include "testing_utilities/is_near.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "numerics/ulp_distance.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_is_near {

using quantities::DebugString;
using quantities::Pow;

template<typename T>
testing::PolymorphicMatcher<IsNearMatcher<ExpectedType<T>>> IsNear(
    T const& expected) {
  return testing::MakePolymorphicMatcher(
      IsNearMatcher<ExpectedType<T>>(expected, /*tolerance=*/1.1));
}

template<typename T>
testing::PolymorphicMatcher<IsNearMatcher<ExpectedType<T>>> IsNear(
    T const& expected,
    double const tolerance) {
  CHECK_LE(1.0, tolerance);
  return testing::MakePolymorphicMatcher(
      IsNearMatcher<ExpectedType<T>>(expected, tolerance));
}

template<typename T>
IsNearMatcher<T>::IsNearMatcher(T const& expected, double const tolerance)
    : expected_(expected),
      low_(std::min(expected * std::sqrt(tolerance),
                    expected / std::sqrt(tolerance))),
      high_(std::max(expected * std::sqrt(tolerance),
                     expected / std::sqrt(tolerance))),
      tolerance_(tolerance) {}

template<typename T>
template<typename Dimensions>
bool IsNearMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  bool const match =  low_ <= actual && actual <= high_;
  if (!match) {
    *listener << "which is not in the range [" << low_ << ", " << high_
              << u8"] and is a factor of √"
              << Pow<2>(std::max(actual / expected_, expected_ / actual))
              << " from the expected value";
  }
  return match;
}

template<typename T>
bool IsNearMatcher<T>::MatchAndExplain(
    double const actual,
    testing::MatchResultListener* listener) const {
  bool const match =  low_ <= actual && actual <= high_;
  if (!match) {
    *listener << "which is not in the range [" << DebugString(low_) << ", "
              << DebugString(high_) << u8"] and is a factor of √"
              << Pow<2>(std::max(actual / expected_, expected_ / actual))
              << " from the expected value";
  }
  return match;
}

template<typename T>
void IsNearMatcher<T>::DescribeTo(std::ostream* out) const {
  *out << "is within ["<< low_
       << ", " << high_ << u8"], i.e., a factor √" << tolerance_
       << " away from " << expected_;
}

template<typename T>
void IsNearMatcher<T>::DescribeNegationTo(std::ostream* out) const {
  *out << "is not within ["<< low_
       << ", " << high_ << u8"], i.e., a factor √" << tolerance_
       << " away from " << expected_;
}

}  // namespace internal_is_near
}  // namespace testing_utilities
}  // namespace principia
