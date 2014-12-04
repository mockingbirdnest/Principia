#pragma once

#include "testing_utilities/componentwise.hpp"

#include <string>

#include "gmock/gmock.h"

using testing::Matcher;

namespace principia {
namespace testing_utilities {

namespace {

template<typename Matcher>
class UnwrapMatcher {
public:
  using matcher = Matcher;
};

template<typename Impl>
class UnwrapMatcher<testing::PolymorphicMatcher<Impl>> {
public:
  using matcher = Impl;
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
  UnwrapMatcher<XMatcher>::matcher(x_matcher_).DescribeTo(out);
  *out << " and ";
  //Matcher<void*>(y_matcher_).DescribeTo(out);
  *out << " and ";
  //Matcher<void*>(z_matcher_).DescribeTo(out);
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
