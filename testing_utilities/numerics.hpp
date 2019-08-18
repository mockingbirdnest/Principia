
#pragma once

#include <cstdint>

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock-matchers.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_numerics {

using quantities::Difference;
using ::testing::Matcher;
using ::testing::MatcherInterface;
using ::testing::MatchResultListener;

template<typename Scalar>
double DoubleValue(Scalar const& scalar);

template<typename T, typename NormType>
NormType AbsoluteError(T const& expected, T const& actual,
                       NormType (T::* norm)() const);

template<typename T, typename NormType, typename NormArg>
NormType AbsoluteError(T const& expected, T const& actual,
                       NormType (*norm)(NormArg const));

// Equivalent to AbsoluteError(expected, actual, &Abs).
double AbsoluteError(double expected, double actual);

// Equivalent to AbsoluteError(expected, actual, &Abs<Dimensions>).
template<typename Dimensions>
quantities::Quantity<Dimensions> AbsoluteError(
    quantities::Quantity<Dimensions> const& expected,
    quantities::Quantity<Dimensions> const& actual);

// Uses |R3Element::Norm|.
template<typename Scalar>
Scalar AbsoluteError(geometry::R3Element<Scalar> const& expected,
                     geometry::R3Element<Scalar> const& actual);

// Uses Multivector.Norm().
template<typename Scalar, typename Frame, int rank>
Scalar AbsoluteError(
    geometry::Multivector<Scalar, Frame, rank> const& expected,
    geometry::Multivector<Scalar, Frame, rank> const& actual);

// Uses the underlying multivector.
template<typename Scalar, typename Frame>
Scalar AbsoluteError(
    geometry::Point<geometry::Multivector<Scalar, Frame, 1>> const& expected,
    geometry::Point<geometry::Multivector<Scalar, Frame, 1>> const& actual);

// Uses the underlying multivector.
template<typename Scalar>
Scalar AbsoluteError(geometry::Point<Scalar> const& expected,
                     geometry::Point<Scalar> const& actual);

template<typename T, typename NormType>
double RelativeError(T const& expected, T const& actual,
                     NormType (T::* norm)() const);

template<typename T, typename NormType, typename NormArg>
double RelativeError(T const& expected, T const& actual,
                     NormType (*norm)(NormArg const));

// Equivalent to RelativeError(expected, actual, &Abs).
double RelativeError(double expected, double actual);

// Equivalent to RelativeError(expected, actual, &Abs<Dimensions>).
template<typename Dimensions>
double RelativeError(quantities::Quantity<Dimensions> const& expected,
                     quantities::Quantity<Dimensions> const& actual);

// Uses |R3Element::Norm|.
template<typename Scalar>
double RelativeError(geometry::R3Element<Scalar> const& expected,
                     geometry::R3Element<Scalar> const& actual);

// Uses Multivector.Norm().
template<typename Scalar, typename Frame, int rank>
double RelativeError(geometry::Multivector<Scalar, Frame, rank> const& expected,
                     geometry::Multivector<Scalar, Frame, rank> const& actual);

template<typename Value>
class DifferenceFromMatcher : public MatcherInterface<Value> {
 public:
  DifferenceFromMatcher(Value const& expected,
                        Matcher<Difference<Value>> const& error_matcher);

  bool MatchAndExplain(Value const& actual,
                       MatchResultListener* listener) const override;
  void DescribeTo(std::ostream* os) const override;
  void DescribeNegationTo(std::ostream* os) const override;

 private:
  Value expected_;
  Matcher<Difference<Value>> error_matcher_;
};

template<typename Value>
class AbsoluteErrorFromMatcher : public MatcherInterface<Value> {
 public:
  using Error =
      decltype(AbsoluteError(std::declval<Value>(), std::declval<Value>()));

  AbsoluteErrorFromMatcher(Value const& expected,
                           Matcher<Error> const& error_matcher);

  bool MatchAndExplain(Value const& actual,
                       MatchResultListener* listener) const override;
  void DescribeTo(std::ostream* os) const override;
  void DescribeNegationTo(std::ostream* os) const override;

 private:
  Value expected_;
  Matcher<Error> error_matcher_;
};

template<typename Value>
class RelativeErrorFromMatcher : public MatcherInterface<Value> {
 public:
  RelativeErrorFromMatcher(Value const& expected,
                           Matcher<double> const& error_matcher);

  bool MatchAndExplain(Value const& actual,
                       MatchResultListener* listener) const override;
  void DescribeTo(std::ostream* os) const override;
  void DescribeNegationTo(std::ostream* os) const override;

 private:
  Value expected_;
  Matcher<double> error_matcher_;
};

template<typename Value>
Matcher<Value> DifferenceFrom(Value const& expected,
                              Matcher<Difference<Value>> const& error_matcher);

template<typename Value, typename ErrorMatcher>
Matcher<Value> AbsoluteErrorFrom(Value const& expected,
                                 ErrorMatcher const& error_matcher);

template<typename Value, typename ErrorMatcher>
Matcher<Value> RelativeErrorFrom(Value const& expected,
                                 Matcher<double> const& error_matcher);

}  // namespace internal_numerics

using internal_numerics::AbsoluteError;
using internal_numerics::AbsoluteErrorFrom;
using internal_numerics::DifferenceFrom;
using internal_numerics::DoubleValue;
using internal_numerics::RelativeError;
using internal_numerics::RelativeErrorFrom;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerics_body.hpp"
