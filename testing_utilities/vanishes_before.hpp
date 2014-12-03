#pragma once

#include <cfloat>
#include <cstdint>

#include <string>

#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
class VanishesBeforeMatcher;

// The matchers below are useful when a computation gives a result whose
// expected value is 0.  Because of cancellations it is unlikely that the
// computed value is exactly zero.  The matchers take a reference value, which
// represents the order of magnitude of the intermediate results that triggered
// the cancellation, and a tolerance expressed as a number of epsilons in the
// error.  More precisely the matcher checks that the actual value is in:
//
//   ]min_epsilons * reference * epsilon, max_epsilons * reference * epsilon]

// The 2-argument version of |VanishesBefore()| should always be preferred as it
// guarantees that the error bound is tight.  The implicit value for
// |min_epsilon| is |max_epsilon| / 2.
template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    double const max_epsilons);

// The 3-argument version of |VanishesBefore()| is exclusively for use when a
// given assertion may have different errors, e.g., because it's in a loop.  It
// doesn't guarantee that the error bound is tight.
template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    double const min_epsilons,
    double const max_epsilons);

template<typename T>
class VanishesBeforeMatcher{
 public:
  explicit VanishesBeforeMatcher(T const& reference,
                                 double const min_epsilons,
                                 double const max_epsilons);
  ~VanishesBeforeMatcher() = default;

  template<typename Dimensions>
  bool MatchAndExplain(quantities::Quantity<Dimensions> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(double const actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  T const reference_;
  double const min_epsilons_;
  double const max_epsilons_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/vanishes_before_body.hpp"
