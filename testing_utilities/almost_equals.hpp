#pragma once

#include <float.h>

#include <string>

#include "gmock/gmock.h"

#include "quantities/dimensionless.hpp"

namespace principia {
namespace test_utilities {

template<typename Scalar>
double DoubleValue(Scalar const& scalar);

double RelativeError(double const expected, double const actual);

template<typename T>
class AlmostEqualsMatcher;

template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(T expected);

template<typename T>
class AlmostEqualsMatcher{
 public:
  explicit AlmostEqualsMatcher(T expected);
  ~AlmostEqualsMatcher() = default;

  template<typename Dimensions>
  virtual bool MatchAndExplain(quantities::Quantity<Dimensions> actual,
                               testing::MatchResultListener * listener) const;
  virtual bool MatchAndExplain(Dimensionless actual,
                               testing::MatchResultListener * listener) const;

  virtual void DescribeTo(std::ostream* os) const;
  virtual void DescribeNegationTo(std::ostream* os) const;

 private:
  T const expected_;
};

template<typename Scalar>
class AlmostVanishesBeforeMatcher : public testing::MatcherInterface<Scalar> {
 public:
  explicit AlmostVanishesBeforeMatcher(Scalar input_magnitude);
  ~AlmostVanishesBeforeMatcher() = default;

  virtual bool MatchAndExplain(Scalar actual,
                               testing::MatchResultListener * listener) const;
  virtual void DescribeTo(std::ostream* os) const;
  virtual void DescribeNegationTo(std::ostream* os) const;

 private:
  Scalar const input_magnitude_;
};

template<typename Scalar>
class ApproximatesMatcher : public testing::MatcherInterface<Scalar> {
 public:
  ApproximatesMatcher(Scalar expected, quantities::Dimensionless relativeError);
  ~ApproximatesMatcher() = default;

  virtual bool MatchAndExplain(Scalar actual,
                               testing::MatchResultListener * listener) const;
  virtual void DescribeTo(std::ostream* os) const;
  virtual void DescribeNegationTo(std::ostream* os) const;

 private:
  Scalar const expected_;
};



MATCHER_P(AlmostVanishesBefore, input_magnitude,
          std::string(negation ? "is not" : "is") + " within " +
          testing::PrintToString(input_magnitude) + " * 4 epsilon of zero.") {
  double const expected_absolute_error =
      DoubleValue(arg_type(input_magnitude)) * 4 * DBL_EPSILON;
  double const actual_absolute_error = std::abs(DoubleValue(arg));
  if (actual_absolute_error <= expected_absolute_error) {
    *result_listener << "the absolute error is " <<
      actual_absolute_error * arg.SIUnit() / (4 * DBL_EPSILON) <<
      " * 4 epsilon";
  }
  return actual_absolute_error <= expected_absolute_error;
}

MATCHER_P2(Approximates, expected, expected_relative_error,
           std::string(negation ? "does not approximate " : "approximates ") +
           testing::PrintToString(expected) + " to within "
           + testing::PrintToString(expected_relative_error)) {
  double const actual_relative_error =
      RelativeError(DoubleValue(arg_type(expected)), DoubleValue(arg));
  if (actual_relative_error <= expected_relative_error) {
    *result_listener << "the relative error is " << actual_relative_error;
  }
  return actual_relative_error <= expected_relative_error;
}


}  // namespace test_utilities
}  // namespace principia

#include "testing_utilities/quantity_matchers_body.hpp"
