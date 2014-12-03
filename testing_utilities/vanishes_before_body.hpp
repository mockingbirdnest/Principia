#pragma once

#include "testing_utilities/numerics.hpp"

#include <float.h>
#include <math.h>
#include <stdint.h>

#include <algorithm>
#include <limits>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"

namespace principia {
namespace testing_utilities {

template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    double const max_epsilons) {
  return testing::MakePolymorphicMatcher(
      VanishesBeforeMatcher<T>(reference, 0.5 * max_epsilons, max_epsilons));
}

template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    double const min_epsilons,
    double const max_epsilons) {
  CHECK_LT(0, min_epsilons);
  CHECK_LT(min_epsilons, max_epsilons);
  return testing::MakePolymorphicMatcher(
      VanishesBeforeMatcher<T>(reference, min_epsilons, max_epsilons));
}

template<typename T>
VanishesBeforeMatcher<T>::VanishesBeforeMatcher(T const& reference,
                                                double const min_epsilons,
                                                double const max_epsilons)
    : reference_(reference),
      min_epsilons_(min_epsilons),
      max_epsilons_(max_epsilons) {}

template<typename T>
template<typename Dimensions>
bool VanishesBeforeMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == reference_) {
    return true;
  }
  double const distance =
      AbsoluteError(DoubleValue(actual), 0.0) /
          (DoubleValue(reference_) * std::numeric_limits<double>::epsilon());
  *listener << "the numbers are separated by " << distance << " epsilons";
  return min_epsilons_ < distance && distance <= max_epsilons_;
}

template<typename T>
bool VanishesBeforeMatcher<T>::MatchAndExplain(
    double const actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == reference_) {
    return true;
  }
  double const distance =
      AbsoluteError(actual, 0.0) /
          (reference_ * std::numeric_limits<double>::epsilon());
  *listener << "the numbers are separated by " << distance << " epsilons";
  return min_epsilons_ < distance && distance <= max_epsilons_;
}

template<typename T>
void VanishesBeforeMatcher<T>::DescribeTo(std::ostream* out) const {
  *out << "is within "<< min_epsilons_
       << " to " << max_epsilons_ << " epsilons of " << reference_;
}

template<typename T>
void VanishesBeforeMatcher<T>::DescribeNegationTo(std::ostream* out) const {
  *out << "is not within " << min_epsilons_
       << " to " << max_epsilons_ << " epsilons of " << reference_;
}

}  // namespace testing_utilities
}  // namespace principia
