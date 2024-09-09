#pragma once

#include <cfloat>
#include <cstdint>

#include <string>

#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace _vanishes_before {
namespace internal {

using namespace principia::quantities::_quantities;

template<typename T>
class VanishesBeforeMatcher;

// The matchers below are useful when a computation gives a result whose
// expected value is zero.  Because of cancellations it is unlikely that the
// computed value is exactly zero.  The matchers take a reference value, which
// represents the order of magnitude of the intermediate results that triggered
// the cancellation, and a tolerance expressed as a number of ulps in the
// error.  More precisely the matcher checks that the `reference` is equal to
// to `actual + reference` to within the specified number of ulps.

// The 2-argument version of `VanishesBefore()` should always be preferred as it
// guarantees that the error bound is tight.
template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    std::int64_t max_ulps);

// The 3-argument version of `VanishesBefore()` is exclusively for use when a
// given assertion may have different errors, e.g., because it's in a loop.  It
// doesn't guarantee that the error bound is tight.
template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    std::int64_t min_ulps,
    std::int64_t max_ulps);

template<typename T>
class VanishesBeforeMatcher final {
 public:
  explicit VanishesBeforeMatcher(T const& reference,
                                 std::int64_t min_ulps,
                                 std::int64_t max_ulps);

  template<typename Dimensions>
  bool MatchAndExplain(Quantity<Dimensions> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(double actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  T const reference_;
  std::int64_t const min_ulps_;
  std::int64_t const max_ulps_;
};

}  // namespace internal

using internal::VanishesBefore;

}  // namespace _vanishes_before
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/vanishes_before_body.hpp"
