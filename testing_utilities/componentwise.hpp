#pragma once

#include <cfloat>
#include <cstdint>

#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

using principia::geometry::R3Element;

namespace principia {
namespace testing_utilities {

template<typename XMatcher, typename YMatcher, typename ZMatcher>
class ComponentwiseMatcher;

template<typename XMatcher, typename YMatcher, typename ZMatcher>
testing::PolymorphicMatcher<ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher);

template<typename XMatcher, typename YMatcher, typename ZMatcher>
class ComponentwiseMatcher {
 public:
  explicit ComponentwiseMatcher(XMatcher const& x_matcher,
                                YMatcher const& y_matcher,
                                ZMatcher const& z_matcher);
  ~ComponentwiseMatcher() = default;

  template<typename Scalar>
  bool MatchAndExplain(geometry::R3Element<Scalar> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Bivector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  XMatcher const& x_matcher_;
  YMatcher const& y_matcher_;
  ZMatcher const& z_matcher_;
};

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/componentwise_body.hpp"
