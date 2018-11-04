
#pragma once

#include "testing_utilities/componentwise.hpp"

#include <string>

#include "base/not_constructible.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_componentwise {

using base::not_constructible;
using geometry::Point;
using quantities::Quantity;
using ::testing::MakeMatcher;
using ::testing::SafeMatcherCast;

template<typename T1Matcher, typename T2Matcher>
testing::PolymorphicMatcher<ComponentwiseMatcher2<T1Matcher, T2Matcher>>
Componentwise(T1Matcher const& t1_matcher,
              T2Matcher const& t2_matcher) {
  return testing::MakePolymorphicMatcher(
      ComponentwiseMatcher2<T1Matcher, T2Matcher>(t1_matcher, t2_matcher));
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
testing::PolymorphicMatcher<ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher) {
  return testing::MakePolymorphicMatcher(
      ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>(
          x_matcher, y_matcher, z_matcher));
}

template<typename T1Matcher, typename T2Matcher>
ComponentwiseMatcher2<T1Matcher, T2Matcher>::ComponentwiseMatcher2(
    T1Matcher const& t1_matcher,
    T2Matcher const& t2_matcher)
    : t1_matcher_(t1_matcher),
      t2_matcher_(t2_matcher) {}

template<typename T1Matcher, typename T2Matcher>
template<typename PairType>
ComponentwiseMatcher2<T1Matcher, T2Matcher>::operator Matcher<PairType>()
    const {
  return MakeMatcher(
      new ComponentwiseMatcher2Impl<PairType>(t1_matcher_, t2_matcher_));
}

template<typename T1Matcher, typename T2Matcher, typename T3Matcher>
ComponentwiseMatcher3<T1Matcher, T2Matcher, T3Matcher>::ComponentwiseMatcher3(
    T1Matcher const& t1_matcher,
    T2Matcher const& t2_matcher,
    T3Matcher const& t3_matcher)
    : t1_matcher_(t1_matcher),
      t2_matcher_(t2_matcher),
      t3_matcher_(t3_matcher) {}

template<typename T1Matcher, typename T2Matcher, typename T3Matcher>
template<typename TripleType>
ComponentwiseMatcher3<T1Matcher, T2Matcher, T3Matcher>::
operator Matcher<TripleType>() const {
  return MakeMatcher(new ComponentwiseMatcher3Impl<TripleType>(
      t1_matcher_, t2_matcher_, t3_matcher_));
}

template<typename T1, typename T2>
template<typename T1Matcher, typename T2Matcher>
ComponentwiseMatcher2Impl<geometry::Pair<T1, T2>>::ComponentwiseMatcher2Impl(
    T1Matcher const& t1_matcher,
    T2Matcher const& t2_matcher)
    : t1_matcher_(SafeMatcherCast<T1>(t1_matcher)),
      t2_matcher_(SafeMatcherCast<T2>(t2_matcher)) {}

template<typename T1, typename T2>
bool ComponentwiseMatcher2Impl<geometry::Pair<T1, T2>>::MatchAndExplain(
    geometry::Pair<T1, T2> const& actual,
    testing::MatchResultListener* listener) const {
  bool const t1_matches = t1_matcher_.MatchAndExplain(actual.t1_, listener);
  if (!t1_matches) {
    *listener << " in the first element; ";
  }
  bool const t2_matches = t2_matcher_.MatchAndExplain(actual.t2_, listener);
  if (!t2_matches) {
    *listener << " in the second element; ";
  }
  return t1_matches && t2_matches;
}

template<typename T1, typename T2>
void ComponentwiseMatcher2Impl<geometry::Pair<T1, T2>>::DescribeTo(
    std::ostream* out) const {
  *out << "first element ";
  t1_matcher_.DescribeTo(out);
  *out << " and second element ";
  t2_matcher_.DescribeTo(out);
}

template<typename T1, typename T2>
void ComponentwiseMatcher2Impl<geometry::Pair<T1, T2>>::DescribeNegationTo(
    std::ostream* out) const {
  *out << "first element ";
  t1_matcher_.DescribeNegationTo(out);
  *out << " or second element ";
  t2_matcher_.DescribeNegationTo(out);
}

template<typename T1Matcher, typename T2Matcher>
template<typename T1, typename T2>
bool ComponentwiseMatcher2<T1Matcher, T2Matcher>::MatchAndExplain(
    geometry::Pair<T1, T2> const& actual,
    testing::MatchResultListener* listener) const {
  bool const t1_matches = Matcher<T1>(t1_matcher_).MatchAndExplain(
                              actual.t1_, listener);
  if (!t1_matches) {
    *listener << " in the t1 coordinate; ";
  }
  bool const t2_matches = Matcher<T2>(t2_matcher_).MatchAndExplain(
                              actual.t2_, listener);
  if (!t2_matches) {
    *listener << " in the t2 coordinate; ";
  }
  return t1_matches && t2_matches;
}

template<typename T1Matcher, typename T2Matcher>
template<typename Scalar, typename Frame>
bool ComponentwiseMatcher2<T1Matcher, T2Matcher>::MatchAndExplain(
    geometry::RP2Point<Scalar, Frame> const& actual,
    testing::MatchResultListener* listener) const {
  bool const x_matches = Matcher<Scalar>(t1_matcher_).MatchAndExplain(
                              actual.x(), listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches = Matcher<Scalar>(t2_matcher_).MatchAndExplain(
                              actual.y(), listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  return x_matches && y_matches;
}

template<typename T1Matcher, typename T2Matcher>
void ComponentwiseMatcher2<T1Matcher, T2Matcher>::DescribeTo(
    std::ostream* out) const {
  *out << "t1 ";
  Matcher<typename MatcherParameterType<T1Matcher>::type>(
      t1_matcher_).DescribeTo(out);
  *out << " and t2 ";
  Matcher<typename MatcherParameterType<T2Matcher>::type>(
      t2_matcher_).DescribeTo(out);
}

template<typename T1Matcher, typename T2Matcher>
void ComponentwiseMatcher2<T1Matcher, T2Matcher>::DescribeNegationTo(
    std::ostream* out) const {
  *out << "t2 ";
  Matcher<typename MatcherParameterType<T1Matcher>::type>(
      t1_matcher_).DescribeNegationTo(out);
  *out << " or t2 ";
  Matcher<typename MatcherParameterType<T2Matcher>::type>(
      t2_matcher_).DescribeNegationTo(out);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
template<typename Scalar>
bool ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>::MatchAndExplain(
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
bool ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>::
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
bool ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>::
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
void ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>::DescribeTo(
    std::ostream* out) const {
  *out << "x ";
  Matcher<typename MatcherParameterType<XMatcher>::type>(
      x_matcher_).DescribeTo(out);
  *out << " and y ";
  Matcher<typename MatcherParameterType<YMatcher>::type>(
      y_matcher_).DescribeTo(out);
  *out << " and z ";
  Matcher<typename MatcherParameterType<ZMatcher>::type>(
      z_matcher_).DescribeTo(out);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
void ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>::DescribeNegationTo(
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

}  // namespace internal_componentwise
}  // namespace testing_utilities
}  // namespace principia
