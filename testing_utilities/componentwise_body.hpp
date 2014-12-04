#pragma once

#include "testing_utilities/componentwise.hpp"

#include <string>

#include "gmock/gmock.h"

using testing::Matcher;

namespace principia {
namespace testing_utilities {

namespace {

template<typename T>
class DescribeHelper {
 public:
  static T const& Cast(T const& t) { return t; }
};

template<typename Impl>
class DescribeHelper<testing::PolymorphicMatcher<Impl>> {
 public:
  static Impl const& Cast(testing::PolymorphicMatcher<Impl> const& m) { return m.impl(); }
};

}  // namespace

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
  return Matcher<Scalar>(x_matcher_).MatchAndExplain(actual.x, listener) &&
         Matcher<Scalar>(y_matcher_).MatchAndExplain(actual.y, listener) &&
         Matcher<Scalar>(z_matcher_).MatchAndExplain(actual.z, listener);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
void ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::DescribeTo(
    std::ostream* out) const {
  DescribeHelper<XMatcher>::Cast(x_matcher_).DescribeTo(out);
  *out << " and ";
  DescribeHelper<YMatcher>::Cast(y_matcher_).DescribeTo(out);
  *out << " and ";
  DescribeHelper<ZMatcher>::Cast(z_matcher_).DescribeTo(out);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
void ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::DescribeNegationTo(
    std::ostream* out) const {
  //Matcher<void*>(x_matcher_).DescribeNegationTo(out);
  *out << " or ";
  //Matcher<void*>(y_matcher_).DescribeNegationTo(out);
  *out << " or ";
  //Matcher<void*>(z_matcher_).DescribeNegationTo(out);
}

}  // namespace testing_utilities
}  // namespace principia
