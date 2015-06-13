#pragma once

#include "testing_utilities/almost_equals.hpp"

#include <float.h>
#include <math.h>
#include <stdint.h>

#include <algorithm>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected,
    std::int64_t const max_ulps) {
  return testing::MakePolymorphicMatcher(
      AlmostEqualsMatcher<T>(expected, max_ulps, max_ulps));
}

template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected,
    std::int64_t const min_ulps,
    std::int64_t const max_ulps) {
  CHECK_LE(min_ulps, max_ulps);
  return testing::MakePolymorphicMatcher(
      AlmostEqualsMatcher<T>(expected, min_ulps, max_ulps));
}

template<typename T>
AlmostEqualsMatcher<T>::AlmostEqualsMatcher(T const& expected,
                                            std::int64_t const min_ulps,
                                            std::int64_t const max_ulps)
    : expected_(expected),
      min_ulps_(min_ulps),
      max_ulps_(max_ulps) {}

template<typename T>
template<typename Dimensions>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    quantities::Quantity<Dimensions> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return MatchAndExplainIdentical(listener);
  }
  std::int64_t const distance = ULPDistance(DoubleValue(actual),
                                            DoubleValue(expected_));
  bool const match =  min_ulps_ <= distance && distance <= max_ulps_;
  if (!match) {
    *listener << "the numbers are separated by " << distance << " ULPs";
  }
  return match;
}

template<typename T>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    double const actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return MatchAndExplainIdentical(listener);
  }
  std::int64_t const distance = ULPDistance(actual, expected_);
  bool const match =  min_ulps_ <= distance && distance <= max_ulps_;
  if (!match) {
    *listener << "the numbers are separated by " << distance << " ULPs";
  }
  return match;
}

template<typename T>
template<typename Scalar>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::R3Element<Scalar> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return MatchAndExplainIdentical(listener);
  }
  std::int64_t const x_distance = ULPDistance(DoubleValue(actual.x),
                                              DoubleValue(expected_.x));
  std::int64_t const y_distance = ULPDistance(DoubleValue(actual.y),
                                              DoubleValue(expected_.y));
  std::int64_t const z_distance = ULPDistance(DoubleValue(actual.z),
                                              DoubleValue(expected_.z));
  std::int64_t const max_distance =
      std::max({x_distance, y_distance, z_distance});
  bool const x_is_max = x_distance == max_distance;
  bool const y_is_max = y_distance == max_distance;
  bool const z_is_max = z_distance == max_distance;
  bool const matches = min_ulps_ <= max_distance && max_distance <= max_ulps_;
  if (!matches) {
    *listener << "the following components are not within "
              << min_ulps_ << " to " << max_ulps_ << " ULPs: "
              << (x_is_max ? "" : "x, ") << (y_is_max ? "" : "y, ")
              << (z_is_max ? "" : "z, ");
    *listener << "the components differ by the following numbers of ULPs: x: "
              << x_distance << ", y: " << y_distance << ", z: " << z_distance;
  }
  return matches;
}

template<typename T>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::Quaternion const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return MatchAndExplainIdentical(listener);
  }
  std::int64_t const w_distance = ULPDistance(
      DoubleValue(actual.real_part()),
      DoubleValue(expected_.real_part()));
  std::int64_t const x_distance = ULPDistance(
      DoubleValue(actual.imaginary_part().x),
      DoubleValue(expected_.imaginary_part().x));
  std::int64_t const y_distance = ULPDistance(
      DoubleValue(actual.imaginary_part().y),
      DoubleValue(expected_.imaginary_part().y));
  std::int64_t const z_distance = ULPDistance(
      DoubleValue(actual.imaginary_part().z),
      DoubleValue(expected_.imaginary_part().z));
  bool const w_matches = w_distance <= max_ulps_;
  bool const x_matches = x_distance <= max_ulps_;
  bool const y_matches = y_distance <= max_ulps_;
  bool const z_matches = z_distance <= max_ulps_;
  bool const matches = w_matches && x_matches && y_matches && z_matches;
  if (!matches) {
    *listener << "the following components differ by more than " << max_ulps_
              << " ULPs: "
              << (w_matches ? "" : "w, ") << (x_matches ? "" : "x, ")
              << (y_matches ? "" : "y, ") << (z_matches ? "" : "z, ");
    *listener << "the components differ by the following numbers of ULPs: w: "
              << w_distance << ", x: " << x_distance
              << ", y: " << y_distance << ", z: " << z_distance;
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
    return MatchAndExplainIdentical(listener);
  }
  return AlmostEqualsMatcher<geometry::R3Element<Scalar>>(
      expected_.coordinates(),
      min_ulps_,
      max_ulps_).MatchAndExplain(actual.coordinates(), listener);
}

template<typename T>
template<typename Scalar, typename Frame>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::Bivector<Scalar, Frame> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return MatchAndExplainIdentical(listener);
  }
  return AlmostEqualsMatcher<geometry::R3Element<Scalar>>(
      expected_.coordinates(),
      min_ulps_,
      max_ulps_).MatchAndExplain(actual.coordinates(), listener);
}

template<typename T>
template<typename Scalar, typename Frame>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::Trivector<Scalar, Frame> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return MatchAndExplainIdentical(listener);
  }
  return AlmostEqualsMatcher<Scalar>(expected_.coordinates(),
                                     min_ulps_,
                                     max_ulps_).MatchAndExplain(
                                         actual.coordinates(),
                                         listener);
}

template<typename T>
template<typename Vector>
bool AlmostEqualsMatcher<T>::MatchAndExplain(
    geometry::Point<Vector> const& actual,
    testing::MatchResultListener* listener) const {
  // Check that the types are equality-comparable up to implicit casts.
  if (actual == expected_) {
    return MatchAndExplainIdentical(listener);
  }
  geometry::Point<Vector> const origin;
  return AlmostEqualsMatcher<Vector>(expected_ - origin,
                                     min_ulps_,
                                     max_ulps_).MatchAndExplain(
                                         actual - origin,
                                         listener);
}

template<typename T>
void AlmostEqualsMatcher<T>::DescribeTo(std::ostream* out) const {
  *out << "is within "<< min_ulps_
       << " to " << max_ulps_ << " ULPs of " << expected_;
}

template<typename T>
void AlmostEqualsMatcher<T>::DescribeNegationTo(std::ostream* out) const {
  *out << "is not within " << min_ulps_
       << " to " << max_ulps_ << " ULPs of " << expected_;
}

template<typename T>
bool AlmostEqualsMatcher<T>::MatchAndExplainIdentical(
    testing::MatchResultListener* listener) const {
  if (min_ulps_ == 0) {
    return true;
  } else {
    *listener << "the numbers are identical";
    return false;
  }
}

}  // namespace testing_utilities
}  // namespace principia
