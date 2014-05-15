#pragma once

#include <cstdint>

#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
testing::PolymorphicMatcher<AlmostVanishesBeforeMatcher<T>>
AlmostVanishesBefore(T const& input_magnitude,
                     std::int64_t const max_ulps = 4) {
  return testing::MakePolymorphicMatcher(
      AlmostVanishesBeforeMatcher(input_magnitude, max_ulps))
}


template<typename T>
AlmostVanishesBeforeMatcher<T>::AlmostVanishesBeforeMatcher(
    T const& expected,
    int64_t const max_ulps)
    : expected_(expected),
      max_ulps_(max_ulps) {}

template<typename T>
template<typename Dimensions>
bool AlmostVanishesBeforeMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  return AlmostEqualsMatcher(
        actual + DoubleValue(input_magnitude) * actual.SIUnit).MatchAndExplain(
            Quantity<Dimensions>(), listener);
}

template<typename T>
bool AlmostVanishesBeforeMatcher<T>::MatchAndExplain(
    quantities::Dimensionless const& actual,
    testing::MatchResultListener* listener) const {
  return AlmostEqualsMatcher(
        actual + DoubleValue(input_magnitude)).MatchAndExplain(0, listener);
}

template<typename T>
template<typename Scalar>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::R3Element<Scalar> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return true;
  }
  int64_t const x_distance =
    testing::internal::Double(DoubleValue(actual.x)).AlmostEquals(
        testing::internal::Double(DoubleValue(expected_.x)));
  int64_t const y_distance =
    testing::internal::Double(DoubleValue(actual.y)).AlmostEquals(
        testing::internal::Double(DoubleValue(expected_.y)));
  int64_t const z_distance =
    testing::internal::Double(DoubleValue(actual.z)).AlmostEquals(
        testing::internal::Double(DoubleValue(expected_.z)));
  bool const x_matches = x_distance <= max_ulps_;
  bool const y_matches = y_distance <= max_ulps_;
  bool const z_matches = z_distance <= max_ulps_;
  bool const matches = x_matches && y_matches && z_matches;
  if (!matches) {
    *listener << "the following components differ by more than " << max_ulps_
              << " ULPs: " << (x_matches ? "" : "x, ") 
              << (y_matches ? "" : "y, ") << (z_matches ? "" : "z, ")
              << "the components differ by the following numbers of ULPs: x: "
              << x_distance << ", y: " << y_distance << ", z: " << z_distance;
  }
  return matches;
}

template<typename T>
template<typename Scalar, typename Frame>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::Vector<Scalar, Frame> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return true;
  }
  return AlmostEqualsMatcher<geometry::R3Element<Scalar>>(
      expected_.coordinates(), 
      max_ulps_).MatchAndExplain(actual.coordinates(), listener);
}

template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeTo(std::ostream* os) const {
  *os << "is within 4 ULPs of " << expected_;
}

template<typename Scalar>
void AlmostEqualsMatcher<Scalar>::DescribeNegationTo(std::ostream* os) const {
  *os << "is not within 4 ULPs of " << expected_;
}


}  // namespace testing_utilities
}  // namespace principia
