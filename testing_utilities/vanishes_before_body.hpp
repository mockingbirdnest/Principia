
#pragma once

#include "testing_utilities/vanishes_before.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

#include "gmock/gmock.h"
#include "numerics/ulp_distance.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_vanishes_before {

using numerics::ULPDistance;

template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    std::int64_t const max_ulps) {
  return testing::MakePolymorphicMatcher(
      VanishesBeforeMatcher<T>(reference, max_ulps, max_ulps));
}

template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    std::int64_t const min_ulps,
    std::int64_t const max_ulps) {
  CHECK_LE(min_ulps, max_ulps);
  return testing::MakePolymorphicMatcher(
      VanishesBeforeMatcher<T>(reference, min_ulps, max_ulps));
}

template<typename T>
VanishesBeforeMatcher<T>::VanishesBeforeMatcher(T const& reference,
                                                std::int64_t const min_ulps,
                                                std::int64_t const max_ulps)
    : reference_(reference),
      min_ulps_(min_ulps),
      max_ulps_(max_ulps) {}

template<typename T>
template<typename Dimensions>
bool VanishesBeforeMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  std::int64_t const distance =
      ULPDistance(DoubleValue(reference_), DoubleValue(actual + reference_));
  bool const matches = min_ulps_ <= distance && distance <= max_ulps_;
  if (!matches) {
    *listener << "the numbers are separated by " << distance << " ulps";
  }
  return matches;
}

template<typename T>
bool VanishesBeforeMatcher<T>::MatchAndExplain(
    double const actual,
    testing::MatchResultListener* listener) const {
  std::int64_t const distance = ULPDistance(reference_, actual + reference_);
  bool const matches = min_ulps_ <= distance && distance <= max_ulps_;
  if (!matches) {
    *listener << "the numbers are separated by " << distance << " ulps";
  }
  return matches;
}

template<typename T>
void VanishesBeforeMatcher<T>::DescribeTo(std::ostream* out) const {
  *out << "vanishes before " << reference_ << " to within "
       << min_ulps_ << " to " << max_ulps_ << " ULPs";
}

template<typename T>
void VanishesBeforeMatcher<T>::DescribeNegationTo(std::ostream* out) const {
  *out << "does not vanish before " << reference_ << " to within "
       << min_ulps_ << " to " << max_ulps_ << " ULP";
}

}  // namespace internal_vanishes_before
}  // namespace testing_utilities
}  // namespace principia
