#pragma once

#include "testing_utilities/componentwise.hpp"

#include <string>

#include "gmock/gmock.h"

namespace principia {
namespace testing_utilities {

template<typename XMatcher, typename YMatcher, typename ZMatcher>
testing::PolymorphicMatcher<ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher) {
  return testing::MakePolymorphicMatcher(
      ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>(
          x_matcher, y_matcher, z_matcher));
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::ComponentwiseMatcher(
    XMatcher const& x_matcher,
    YMatcher const& y_matcher,
    ZMatcher const& z_matcher)
    : x_matcher_(x_matcher),
      y_matcher_(y_matcher),
      z_matcher_(z_matcher) {}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
template<typename Scalar>
bool ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::MatchAndExplain(
    geometry::R3Element<Scalar> const& actual,
    testing::MatchResultListener* listener) const {
  return x_matcher_.MatchAndExplain(actual.x, listener) &&
         y_matcher_.MatchAndExplain(actual.y, listener) &&
         z_matcher_.MatchAndExplain(actual.z, listener);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
void ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::DescribeTo(
    std::ostream* out) const {
  x_matcher_.DescribeTo(out);
  *out << " and ";
  y_matcher_.DescribeTo(out);
  *out << " and ";
  z_matcher_.DescribeTo(out);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
void ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::DescribeNegationTo(
    std::ostream* out) const {
  x_matcher_.DescribeNegationTo(out);
  *out << " or ";
  y_matcher_.DescribeNegationTo(out);
  *out << " or ";
  z_matcher_.DescribeNegationTo(out);
}

}  // namespace testing_utilities
}  // namespace principia
