#pragma once

#include "testing_utilities/numerics_matchers.hpp"

#include <ostream>
#include <utility>

#include "testing_utilities/numerics.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_numerics_matchers {

using ::testing::MakeMatcher;
using ::testing::MatcherInterface;
using ::testing::MatchResultListener;

template<typename Value>
class DifferenceFromMatcher : public MatcherInterface<Value> {
 public:
  DifferenceFromMatcher(Value const& expected,
                        Matcher<Difference<Value>> error_matcher);

  bool MatchAndExplain(Value actual,
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
                           Matcher<Error> error_matcher);

  bool MatchAndExplain(Value actual,
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
                           Matcher<double> error_matcher);

  bool MatchAndExplain(Value actual,
                       MatchResultListener* listener) const override;
  void DescribeTo(std::ostream* os) const override;
  void DescribeNegationTo(std::ostream* os) const override;

 private:
  Value expected_;
  Matcher<double> error_matcher_;
};

template<typename Value>
DifferenceFromMatcher<Value>::DifferenceFromMatcher(
    Value const& expected,
    Matcher<Difference<Value>> error_matcher)
    : expected_(expected),
      error_matcher_(std::move(error_matcher)) {}

template<typename Value>
bool DifferenceFromMatcher<Value>::MatchAndExplain(
    Value const actual,
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
    Matcher<Error> error_matcher)
    : expected_(expected),
      error_matcher_(std::move(error_matcher)) {}

template<typename Value>
bool AbsoluteErrorFromMatcher<Value>::MatchAndExplain(
    Value const actual,
    MatchResultListener* listener) const {
  Error const error = AbsoluteError(expected_, actual);
  *listener << "whose absolute error from the expected value is " << error
            << " ";
  return error_matcher_.MatchAndExplain(error, listener);
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
    Matcher<double> error_matcher)
    : expected_(expected),
      error_matcher_(std::move(error_matcher)) {}

template<typename Value>
bool RelativeErrorFromMatcher<Value>::MatchAndExplain(
    Value const actual,
    MatchResultListener* listener) const {
  double const error = RelativeError(expected_, actual);
  *listener << "whose relative error from the expected value is " << error
            << " ";
  return error_matcher_.MatchAndExplain(error, listener);
}

template<typename Value>
void RelativeErrorFromMatcher<Value>::DescribeTo(std::ostream* os) const {
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
Matcher<Value> DifferenceFrom(
    Value const& expected,
    Matcher<Difference<Value>> const& error_matcher) {
  return MakeMatcher(new DifferenceFromMatcher<Value>(expected, error_matcher));
}

template<typename Value, typename ErrorMatcher>
Matcher<Value> AbsoluteErrorFrom(Value const& expected,
                                 ErrorMatcher const& error_matcher) {
  return MakeMatcher(
      new AbsoluteErrorFromMatcher<Value>(expected, error_matcher));
}

template<typename Value>
Matcher<Value> RelativeErrorFrom(Value const& expected,
                                 Matcher<double> const& error_matcher) {
  return MakeMatcher(
      new RelativeErrorFromMatcher<Value>(expected, error_matcher));
}

}  // namespace internal_numerics_matchers
}  // namespace testing_utilities
}  // namespace principia
