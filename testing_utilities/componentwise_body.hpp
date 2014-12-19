#pragma once

#include "testing_utilities/componentwise.hpp"

#include <string>

#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

using principia::quantities::Quantity;
using testing::Matcher;

namespace principia {
namespace testing_utilities {

namespace {

// In order to call |Describe...To| on the various matchers we need to convert
// them to some |Matcher<T>|.  However, it is quite difficult to figure out
// what |T| is because the class |ComponentwiseMatcher| and the factory
// |Componentwise| cannot take it as a parameter (the template deduction would
// fail and the usages would get very ugly).  To make things worse, the actual
// type of the matcher depends on whether it polymorphic, monomorphic or
// some other internal helper.  We obtain |T| by peeling away the layers of
// templates around it.

template<typename T>
class MatcherParameterType {
 public:
  using type = T;
};

template<typename T, template<typename> class U>
class MatcherParameterType<U<T>> {
 public:
  using type = typename MatcherParameterType<T>::type;
};

// And now a case that we *don't* want to peel away.  Yes, this smells a bit.
template<typename T>
class MatcherParameterType<Quantity<T>> {
 public:
  using type = Quantity<T>;
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
  bool const x_matches =  Matcher<Scalar>(x_matcher_).MatchAndExplain(
                              actual.x, listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches =  Matcher<Scalar>(y_matcher_).MatchAndExplain(
                              actual.y, listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  bool const z_matches =  Matcher<Scalar>(z_matcher_).MatchAndExplain(
                              actual.z, listener);
  if (!z_matches) {
    *listener << " in the z coordinate; ";
  }
  return x_matches && y_matches && z_matches;
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
template<typename Scalar, typename Frame>
bool ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::
MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                testing::MatchResultListener* listener) const {
  bool const x_matches =  Matcher<Scalar>(x_matcher_).MatchAndExplain(
                              actual.coordinates().x, listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches =  Matcher<Scalar>(y_matcher_).MatchAndExplain(
                              actual.coordinates().y, listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  bool const z_matches =  Matcher<Scalar>(z_matcher_).MatchAndExplain(
                              actual.coordinates().z, listener);
  if (!z_matches) {
    *listener << " in the z coordinate; ";
  }
  return x_matches && y_matches && z_matches;
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
template<typename Scalar, typename Frame>
bool ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::
MatchAndExplain(geometry::Bivector<Scalar, Frame> const& actual,
                testing::MatchResultListener* listener) const {
  bool const x_matches =  Matcher<Scalar>(x_matcher_).MatchAndExplain(
                              actual.coordinates().x, listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches =  Matcher<Scalar>(y_matcher_).MatchAndExplain(
                              actual.coordinates().y, listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  bool const z_matches =  Matcher<Scalar>(z_matcher_).MatchAndExplain(
                              actual.coordinates().z, listener);
  if (!z_matches) {
    *listener << " in the z coordinate; ";
  }
  return x_matches && y_matches && z_matches;
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
void ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::DescribeTo(
    std::ostream* out) const {
  *out << "x ";
  Matcher<typename MatcherParameterType<XMatcher>::type>(x_matcher_).DescribeTo(out);
  *out << " and y ";
  Matcher<typename MatcherParameterType<YMatcher>::type>(y_matcher_).DescribeTo(out);
  *out << " and z ";
  Matcher<typename MatcherParameterType<ZMatcher>::type>(z_matcher_).DescribeTo(out);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
void ComponentwiseMatcher<XMatcher, YMatcher, ZMatcher>::DescribeNegationTo(
    std::ostream* out) const {
  *out << "x ";
  Matcher<typename MatcherParameterType<XMatcher>::type>(
      x_matcher_).DescribeNegationTo(out);
  *out << " or y ";
  Matcher<typename MatcherParameterType<YMatcher>::type>(
      y_matcher_).DescribeNegationTo(out);
  *out << " or z ";
  Matcher<typename MatcherParameterType<ZMatcher>::type>(
      z_matcher_).DescribeNegationTo(out);
}

}  // namespace testing_utilities
}  // namespace principia
