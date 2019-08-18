
#pragma once

#include "testing_utilities/numerics.hpp"

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>

namespace principia {
namespace testing_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar) {
  return scalar / quantities::SIUnit<Scalar>();
}

template<typename T, typename NormType>
NormType AbsoluteError(T const& expected, T const& actual,
                       NormType (T::* norm)() const) {
  return ((expected - actual).*norm)();
}

template<typename T, typename NormType, typename NormArg>
NormType AbsoluteError(T const& expected, T const& actual,
                       NormType (*norm)(NormArg const)) {
  return norm(expected - actual);
}

inline double AbsoluteError(double const expected, double const actual) {
  return AbsoluteError(expected, actual, &quantities::Abs<double>);
}

template<typename Dimensions>
quantities::Quantity<Dimensions> AbsoluteError(
    quantities::Quantity<Dimensions> const& expected,
    quantities::Quantity<Dimensions> const& actual) {
  return AbsoluteError(
      expected, actual, &quantities::Abs<quantities::Quantity<Dimensions>>);
}

template<typename Scalar>
Scalar AbsoluteError(geometry::R3Element<Scalar> const& expected,
                     geometry::R3Element<Scalar> const& actual) {
  return AbsoluteError(expected, actual, &geometry::R3Element<Scalar>::Norm);
}

template<typename Scalar, typename Frame, int rank>
Scalar AbsoluteError(
    geometry::Multivector<Scalar, Frame, rank> const& expected,
    geometry::Multivector<Scalar, Frame, rank> const& actual) {
  return AbsoluteError<
      geometry::Multivector<Scalar, Frame, rank>, Scalar>(
          expected,
          actual,
          &geometry::Multivector<Scalar, Frame, rank>::Norm);
}

template<typename Scalar, typename Frame>
Scalar AbsoluteError(
    geometry::Point<geometry::Multivector<Scalar, Frame, 1>> const& expected,
    geometry::Point<geometry::Multivector<Scalar, Frame, 1>> const& actual) {
  geometry::Point<geometry::Multivector<Scalar, Frame, 1>> const origin;
  return AbsoluteError(expected - origin, actual - origin);
}

template<typename Scalar>
Scalar AbsoluteError(geometry::Point<Scalar> const& expected,
                     geometry::Point<Scalar> const& actual) {
  geometry::Point<Scalar> const origin;
  return AbsoluteError(expected - origin, actual - origin);
}

template<typename T, typename NormType>
double RelativeError(T const& expected, T const& actual,
                     NormType (T::* norm)() const) {
  if (expected == actual) {
    return 0;
  } else {
    return ((expected - actual).*norm)() / (expected.*norm)();
  }
}

template<typename T, typename NormType, typename NormArg>
double RelativeError(T const& expected, T const& actual,
                     NormType (*norm)(NormArg const)) {
  if (expected == actual) {
    return 0;
  } else {
    return norm(expected - actual) / norm(expected);
  }
}

inline double RelativeError(double const expected, double const actual) {
  return RelativeError(expected, actual, &quantities::Abs<double>);
}

template<typename Dimensions>
double RelativeError(quantities::Quantity<Dimensions> const& expected,
                     quantities::Quantity<Dimensions> const& actual) {
  return RelativeError(
      expected, actual, &quantities::Abs<quantities::Quantity<Dimensions>>);
}

template<typename Scalar>
double RelativeError(geometry::R3Element<Scalar> const& expected,
                     geometry::R3Element<Scalar> const& actual) {
  return RelativeError(expected, actual, &geometry::R3Element<Scalar>::Norm);
}

template<typename Scalar, typename Frame, int rank>
double RelativeError(geometry::Multivector<Scalar, Frame, rank> const& expected,
                     geometry::Multivector<Scalar, Frame, rank> const& actual) {
  return RelativeError(expected, actual,
                       &geometry::Multivector<Scalar, Frame, rank>::Norm);
}

}  // namespace testing_utilities
}  // namespace principia
