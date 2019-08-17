
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
testing::PolymorphicMatcher<IsNearMatcher<T>> IsNear(
    ApproximateQuantity<T> const& expected) {
  return testing::MakePolymorphicMatcher(IsNearMatcher<T>(expected));
}

template<typename T>
IsNearMatcher<T>::IsNearMatcher(ApproximateQuantity<T> const& expected)
    : expected_(expected) {}

template<typename T>
template<typename Dimensions>
bool IsNearMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  bool const match =  expected_.min() <= actual && actual <= expected_.max();
  if (expected_.has_trivial_unit()) {
    *listener << "which ";
  } else {
    *listener << "which is " << actual / expected_.unit() << " * "
              << expected_.unit() << " and ";
  }
  if (match) {
    *listener << "is near " << expected_ << " (being "
              << expected_.UlpDistance(actual) << " ulps away)";
  } else {
    *listener << "is not near " << expected_ << " (being "
              << expected_.UlpDistance(actual) << " ulps away)";
  }
  return match;
}

template<typename T>
bool IsNearMatcher<T>::MatchAndExplain(
    double const actual,
    testing::MatchResultListener* listener) const {
  bool const match =  expected_.min() <= actual && actual <= expected_.max();
  if (match) {
    *listener << "which is near " << expected_ << " (being "
              << expected_.UlpDistance(actual) << " ulps away)";
  } else {
    *listener << "which is not near " << expected_ << " (being "
              << expected_.UlpDistance(actual) << " ulps away)";
  }
  return match;
}

template<typename T>
void IsNearMatcher<T>::DescribeTo(std::ostream* out) const {
  *out << "is near " << expected_;
}

template<typename T>
void IsNearMatcher<T>::DescribeNegationTo(std::ostream* out) const {
  *out << "is not near " << expected_;
}

}  // namespace internal_is_near
}  // namespace testing_utilities
}  // namespace principia
