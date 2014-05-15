#pragma once

#include <cstdint>

#include "gmock/gmock.h"

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "quantities/dimensionless.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
class AlmostVanishesBeforeMatcher;

template<typename T>
testing::PolymorphicMatcher<AlmostVanishesBeforeMatcher<T>>
AlmostVanishesBefore(T const& input_magnitude,
                     std::int64_t const max_ulps = 4);

template<typename T>
class AlmostVanishesBeforeMatcher{
 public:
  AlmostVanishesBeforeMatcher(T const& input_magnitude,
                              std::int64_t const max_ulps);
  ~AlmostVanishesBeforeMatcher() = default;

  template<typename Dimensions>
  bool MatchAndExplain(quantities::Quantity<Dimensions> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(quantities::Dimensionless const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(geometry::R3Element<Scalar> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* os) const;
  void DescribeNegationTo(std::ostream* os) const;

 private:
  T const input_magnitude_;
  int64_t const max_ulps_;
};

}  // namespace testing_utilities
}  // namespace principia
