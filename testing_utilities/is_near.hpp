#pragma once

#include <string>
#include <type_traits>

#include "gmock/gmock.h"
#include "quantities/quantities.hpp"
#include "testing_utilities/approximate_quantity.hpp"

namespace principia {
namespace testing_utilities {
namespace _is_near {
namespace internal {

using namespace principia::quantities::_quantities;

template<typename T>
class IsNearMatcher;

// Checks that actual is in the range defined by expected.
template<typename T>
testing::PolymorphicMatcher<IsNearMatcher<T>> IsNear(
    ApproximateQuantity<T> const& expected);

template<typename T>
class IsNearMatcher final {
 public:
  explicit IsNearMatcher(ApproximateQuantity<T> expected);

  template<typename Dimensions>
  bool MatchAndExplain(Quantity<Dimensions> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(double actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  ApproximateQuantity<T> const expected_;
};

}  // namespace internal

using internal::IsNear;

}  // namespace _is_near
}  // namespace testing_utilities
}  // namespace principia

namespace principia::testing_utilities {
using namespace principia::testing_utilities::_is_near;
}  // namespace principia::testing_utilities

#include "testing_utilities/is_near_body.hpp"
