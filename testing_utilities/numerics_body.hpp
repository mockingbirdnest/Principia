#pragma once

#include "testing_utilities/numerics.hpp"

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace testing_utilities {
namespace _numerics {
namespace internal {

using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

template<typename Scalar>
double DoubleValue(Scalar const& scalar) {
  return scalar / si::Unit<Scalar>;
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
  return AbsoluteError(expected, actual, &Abs<double>);
}

template<typename Dimensions>
Quantity<Dimensions> AbsoluteError(
    Quantity<Dimensions> const& expected,
    Quantity<Dimensions> const& actual) {
  return AbsoluteError(
      expected, actual, &Abs<Quantity<Dimensions>>);
}

template<typename Scalar>
Scalar AbsoluteError(R3Element<Scalar> const& expected,
                     R3Element<Scalar> const& actual) {
  return AbsoluteError(expected, actual, &R3Element<Scalar>::Norm);
}

template<typename Scalar, typename Frame, int rank>
Scalar AbsoluteError(
    Multivector<Scalar, Frame, rank> const& expected,
    Multivector<Scalar, Frame, rank> const& actual) {
  return AbsoluteError<
      Multivector<Scalar, Frame, rank>, Scalar>(
          expected,
          actual,
          &Multivector<Scalar, Frame, rank>::Norm);
}

template<typename Scalar, typename Frame>
Scalar AbsoluteError(
    Point<Multivector<Scalar, Frame, 1>> const& expected,
    Point<Multivector<Scalar, Frame, 1>> const& actual) {
  Point<Multivector<Scalar, Frame, 1>> const origin;
  return AbsoluteError(expected - origin, actual - origin);
}

template<typename Scalar>
Scalar AbsoluteError(Point<Scalar> const& expected,
                     Point<Scalar> const& actual) {
  Point<Scalar> const origin;
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
  return RelativeError(expected, actual, &Abs<double>);
}

template<typename Dimensions>
double RelativeError(Quantity<Dimensions> const& expected,
                     Quantity<Dimensions> const& actual) {
  return RelativeError(
      expected, actual, &Abs<Quantity<Dimensions>>);
}

template<typename Scalar>
double RelativeError(R3Element<Scalar> const& expected,
                     R3Element<Scalar> const& actual) {
  return RelativeError(expected, actual, &R3Element<Scalar>::Norm);
}

template<typename Scalar, typename Frame, int rank>
double RelativeError(Multivector<Scalar, Frame, rank> const& expected,
                     Multivector<Scalar, Frame, rank> const& actual) {
  return RelativeError(expected, actual,
                       &Multivector<Scalar, Frame, rank>::Norm);
}

}  // namespace internal
}  // namespace _numerics
}  // namespace testing_utilities
}  // namespace principia
