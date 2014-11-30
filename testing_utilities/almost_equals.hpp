#pragma once

#include <cfloat>
#include <cstdint>

#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"

namespace principia {
namespace testing_utilities {

template<typename T>
class AlmostEqualsMatcher;

// The 2-argument version of |AlmostEquals()| should always be preferred as it
// guarantees that the error bound is tight.
template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected,
    std::int64_t const max_ulps);

// The 3-argument version of |AlmostEquals()| is exclusively for use when a
// given assertion may have different errors, e.g., because it's in a loop.  It
// doesn't guarantee that the error bound is tight.  For vectors, it applies
// only to the component with the largest error.
template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected,
    std::int64_t const min_ulps,
    std::int64_t const max_ulps);

template<typename T>
class AlmostEqualsMatcher{
 public:
  explicit AlmostEqualsMatcher(T const& expected,
                               std::int64_t const min_ulps,
                               std::int64_t const max_ulps);
  ~AlmostEqualsMatcher() = default;

  template<typename Dimensions>
  bool MatchAndExplain(quantities::Quantity<Dimensions> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(double const actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(geometry::R3Element<Scalar> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Bivector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Trivector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  T const expected_;
  std::int64_t const min_ulps_;
  std::int64_t const max_ulps_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/almost_equals_body.hpp"
