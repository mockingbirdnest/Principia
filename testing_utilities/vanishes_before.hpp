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

// The 2-argument version of |VanishesBefore()| should always be preferred as it
// guarantees that the error bound is tight.
template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    std::int64_t const max_epsilons);

// The 3-argument version of |VanishesBefore()| is exclusively for use when a
// given assertion may have different errors, e.g., because it's in a loop.  It
// doesn't guarantee that the error bound is tight.
template<typename T>
testing::PolymorphicMatcher<VanishesBeforeMatcher<T>> VanishesBefore(
    T const& reference,
    std::int64_t const min_epsilons,
    std::int64_t const max_epsilons);

template<typename T>
class VanishesBeforeMatcher{
 public:
  explicit VanishesBeforeMatcher(T const& reference,
                                 std::int64_t const min_epsilons,
                                 std::int64_t const max_epsilons);
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
  std::int64_t const min_epsilons_;
  std::int64_t const max_epsilons_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/vanishes_before_body.hpp"
