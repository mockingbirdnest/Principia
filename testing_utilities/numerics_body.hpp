
#pragma once

#include "testing_utilities/numerics.hpp"

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>

namespace principia {
namespace testing_utilities {
namespace internal_numerics {

using ::testing::MakeMatcher;

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

template<typename Value>
DifferenceFromMatcher<Value>::DifferenceFromMatcher(
    Value const& expected,
    Matcher<Difference<Value>> const& error_matcher)
    : expected_(expected), error_matcher_(error_matcher) {}

template<typename Value>
bool DifferenceFromMatcher<Value>::MatchAndExplain(
    Value const& actual,
    MatchResultListener* listener) const {
  Difference<Value> const difference = actual - expected_;
  *listener << "whose difference from the expected value is " << difference
            << " ";
  return error_matcher_.MatchAndExplain(difference, listener);
}

template<typename Value>
void DifferenceFromMatcher<Value>::DescribeTo(std::ostream* os) const {
  *os << "differs from " << expected_ << " by a value that ";
  error_matcher_.DescribeTo(os);
}

template<typename Value>
void DifferenceFromMatcher<Value>::DescribeNegationTo(std::ostream* os) const {
  *os << "differs from " << expected_ << " by a value that ";
  error_matcher_.DescribeNegationTo(os);
}

template<typename Value>
AbsoluteErrorFromMatcher<Value>::AbsoluteErrorFromMatcher(
    Value const& expected,
    Matcher<Error> const& error_matcher)
    : expected_(expected), error_matcher_(error_matcher) {}

template<typename Value>
bool AbsoluteErrorFromMatcher<Value>::MatchAndExplain(
    Value const& actual,
    MatchResultListener* listener) const {
  Error const error = AbsoluteError(expected_, actual);
  *listener << "whose absolute error from the expected value is " << error
            << " ";
  return error_matcher_.MatchAndExplain(difference, listener);
}

template<typename Value>
void AbsoluteErrorFromMatcher<Value>::DescribeTo(std::ostream* os) const {
  *os << "has an absolute error from " << expected_ << " that ";
  error_matcher_.DescribeTo(os);
}

template<typename Value>
void AbsoluteErrorFromMatcher<Value>::DescribeNegationTo(
    std::ostream* os) const {
  *os << "has an absolute error from " << expected_ << " that ";
  error_matcher_.DescribeNegationTo(os);
}

template<typename Value>
RelativeErrorFromMatcher<Value>::RelativeErrorFromMatcher(
    Value const& expected,
    Matcher<double> const& error_matcher)
    : expected_(expected), error_matcher_(error_matcher) {}

template<typename Value>
bool RelativeErrorFromMatcher<Value>::MatchAndExplain(
    Value const& actual,
    MatchResultListener* listener) const {
  double const error = RelativeError(expected_, actual);
  *listener << "whose relative error from the expected value is " << error
            << " ";
  return error_matcher_.MatchAndExplain(difference, listener);
}

template<typename Value>
void RelativeErrorFromMatcher<Value>::DescribeTo(
    std::ostream* os) const {
  *os << "has a relative error from " << expected_ << " that ";
  error_matcher_.DescribeTo(os);
}

template<typename Value>
inline void RelativeErrorFromMatcher<Value>::DescribeNegationTo(
    std::ostream* os) const {
  *os << "has a relative error from " << expected_ << " that ";
  error_matcher_.DescribeNegationTo(os);
}

template<typename Value>
Matcher<Value> DifferenceFrom(Value const& expected,
                              Matcher<Difference<Value>> const& error_matcher) {
  return MakeMatcher(new DifferenceFromMatcher<Value>(expected, error_matcher));
}

template<typename Value, typename ErrorMatcher>
Matcher<Value> AbsoluteErrorFrom(Value const& expected,
                                 ErrorMatcher const& error_matcher) {
  return MakeMatcher(
      new AbsoluteErrorFromMatcher<Value>(expected, error_matcher));
}

template<typename Value, typename ErrorMatcher>
Matcher<Value> RelativeErrorFrom(Value const& expected,
                                 Matcher<double> const& error_matcher) {
  return MakeMatcher(
      new RelativeErrorFromMatcher<Value>(expected, error_matcher));
}

}  // namespace internal_numerics
}  // namespace testing_utilities
}  // namespace principia
