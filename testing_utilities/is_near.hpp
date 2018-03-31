
#pragma once

#include <string>
#include <type_traits>

#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_is_near {

template<typename T>
class IsNearMatcher;

template<typename T>
using ExpectedType = std::conditional_t<std::is_arithmetic<T>::value, double, T>;

// Calls the next function with |tolerance| set to 1.1.
template<typename T>
testing::PolymorphicMatcher<IsNearMatcher<ExpectedType<T>>> IsNear(
    T const& expected);

// Checks that |expected| is in the range
// [expected / √tolerance, expected √tolerance].
template<typename T>
testing::PolymorphicMatcher<IsNearMatcher<ExpectedType<T>>> IsNear(
    T const& expected,
    double tolerance);

template<typename T>
class IsNearMatcher final {
 public:
  IsNearMatcher(T const& expected,
                double tolerance);

  template<typename Dimensions>
  bool MatchAndExplain(quantities::Quantity<Dimensions> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(double actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  T const expected_;
  T const low_;
  T const high_;
  double tolerance_;
};

}  // namespace internal_is_near

using internal_is_near::IsNear;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/is_near_body.hpp"
