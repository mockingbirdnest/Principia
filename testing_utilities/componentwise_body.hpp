
#pragma once

#include "testing_utilities/componentwise.hpp"

#include <string>

#include "gmock/gmock.h"

namespace principia {
namespace testing_utilities {
namespace internal_componentwise {

using ::testing::SafeMatcherCast;

template<typename T1Matcher, typename T2Matcher>
ComponentwiseMatcher2<T1Matcher, T2Matcher>
Componentwise(T1Matcher const& t1_matcher,
              T2Matcher const& t2_matcher) {
  return ComponentwiseMatcher2<T1Matcher, T2Matcher>(t1_matcher, t2_matcher);
}

template<typename XMatcher, typename YMatcher, typename ZMatcher>
ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher) {
  return ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>(
             x_matcher, y_matcher, z_matcher);
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
  return Matcher<PairType>(
      new ComponentwiseMatcher2Impl<PairType const&>(t1_matcher_, t2_matcher_));
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
  return Matcher<TripleType>(new ComponentwiseMatcher3Impl<TripleType const&>(
      t1_matcher_, t2_matcher_, t3_matcher_));
}

template<typename T1, typename T2>
template<typename T1Matcher, typename T2Matcher>
ComponentwiseMatcher2Impl<geometry::Pair<T1, T2> const&>::
ComponentwiseMatcher2Impl(T1Matcher const& t1_matcher,
                          T2Matcher const& t2_matcher)
    : t1_matcher_(SafeMatcherCast<T1>(t1_matcher)),
      t2_matcher_(SafeMatcherCast<T2>(t2_matcher)) {}

template<typename T1, typename T2>
bool ComponentwiseMatcher2Impl<geometry::Pair<T1, T2> const&>::MatchAndExplain(
    geometry::Pair<T1, T2> const& actual,
    MatchResultListener* listener) const {
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
void ComponentwiseMatcher2Impl<geometry::Pair<T1, T2> const&>::DescribeTo(
    std::ostream* out) const {
  *out << "first element ";
  t1_matcher_.DescribeTo(out);
  *out << " and second element ";
  t2_matcher_.DescribeTo(out);
}

template<typename T1, typename T2>
void ComponentwiseMatcher2Impl<geometry::Pair<T1, T2> const&>::
DescribeNegationTo(std::ostream* out) const {
  *out << "first element ";
  t1_matcher_.DescribeNegationTo(out);
  *out << " or second element ";
  t2_matcher_.DescribeNegationTo(out);
}

template<typename Scalar, typename Frame>
template<typename XMatcher, typename YMatcher>
ComponentwiseMatcher2Impl<geometry::RP2Point<Scalar, Frame> const&>::
ComponentwiseMatcher2Impl(XMatcher const& x_matcher,
                          YMatcher const& y_matcher)
    : x_matcher_(SafeMatcherCast<Scalar>(x_matcher)),
      y_matcher_(SafeMatcherCast<Scalar>(y_matcher)) {}

template<typename Scalar, typename Frame>
bool ComponentwiseMatcher2Impl<geometry::RP2Point<Scalar, Frame> const&>::
MatchAndExplain(geometry::RP2Point<Scalar, Frame> const& actual,
                MatchResultListener* listener) const {
  bool const x_matches = x_matcher_.MatchAndExplain(actual.x(), listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches = y_matcher_.MatchAndExplain(actual.y(), listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  return x_matches && y_matches;
}

template<typename Scalar, typename Frame>
void ComponentwiseMatcher2Impl<geometry::RP2Point<Scalar, Frame> const&>::
DescribeTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeTo(out);
  *out << " and y ";
  y_matcher_.DescribeTo(out);
}

template<typename Scalar, typename Frame>
void ComponentwiseMatcher2Impl<geometry::RP2Point<Scalar, Frame> const&>::
DescribeNegationTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeNegationTo(out);
  *out << " or y ";
  y_matcher_.DescribeNegationTo(out);
}

template<typename Scalar>
template<typename XMatcher, typename YMatcher, typename ZMatcher>
ComponentwiseMatcher3Impl<geometry::R3Element<Scalar> const&>::
ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                          YMatcher const& y_matcher,
                          ZMatcher const& z_matcher)
    : x_matcher_(SafeMatcherCast<Scalar>(x_matcher)),
      y_matcher_(SafeMatcherCast<Scalar>(y_matcher)),
      z_matcher_(SafeMatcherCast<Scalar>(z_matcher)) {}

template<typename Scalar>
bool ComponentwiseMatcher3Impl<geometry::R3Element<Scalar> const&>::
MatchAndExplain(geometry::R3Element<Scalar> const& actual,
                MatchResultListener* listener) const {
  bool const x_matches = x_matcher_.MatchAndExplain(actual.x, listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches = y_matcher_.MatchAndExplain(actual.y, listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  bool const z_matches = z_matcher_.MatchAndExplain(actual.z, listener);
  if (!z_matches) {
    *listener << " in the z coordinate; ";
  }
  return x_matches && y_matches && z_matches;
}

template<typename Scalar>
void ComponentwiseMatcher3Impl<geometry::R3Element<Scalar> const&>::
DescribeTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeTo(out);
  *out << " and y ";
  y_matcher_.DescribeTo(out);
  *out << " and z ";
  z_matcher_.DescribeTo(out);
}

template<typename Scalar>
void ComponentwiseMatcher3Impl<geometry::R3Element<Scalar> const&>::
DescribeNegationTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeNegationTo(out);
  *out << " or y ";
  y_matcher_.DescribeNegationTo(out);
  *out << " or z ";
  z_matcher_.DescribeNegationTo(out);
}

template<typename Scalar, typename Frame>
template<typename XMatcher, typename YMatcher, typename ZMatcher>
ComponentwiseMatcher3Impl<geometry::Vector<Scalar, Frame> const&>::
ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                          YMatcher const& y_matcher,
                          ZMatcher const& z_matcher)
    : x_matcher_(SafeMatcherCast<Scalar>(x_matcher)),
      y_matcher_(SafeMatcherCast<Scalar>(y_matcher)),
      z_matcher_(SafeMatcherCast<Scalar>(z_matcher)) {}

template<typename Scalar, typename Frame>
bool ComponentwiseMatcher3Impl<geometry::Vector<Scalar, Frame> const&>::
MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                MatchResultListener* listener) const {
  bool const x_matches =
      x_matcher_.MatchAndExplain(actual.coordinates().x, listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches =
      y_matcher_.MatchAndExplain(actual.coordinates().y, listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  bool const z_matches =
      z_matcher_.MatchAndExplain(actual.coordinates().z, listener);
  if (!z_matches) {
    *listener << " in the z coordinate; ";
  }
  return x_matches && y_matches && z_matches;
}

template<typename Scalar, typename Frame>
void ComponentwiseMatcher3Impl<geometry::Vector<Scalar, Frame> const&>::
DescribeTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeTo(out);
  *out << " and y ";
  y_matcher_.DescribeTo(out);
  *out << " and z ";
  z_matcher_.DescribeTo(out);
}

template<typename Scalar, typename Frame>
void ComponentwiseMatcher3Impl<geometry::Vector<Scalar, Frame> const&>::
DescribeNegationTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeNegationTo(out);
  *out << " or y ";
  y_matcher_.DescribeNegationTo(out);
  *out << " or z ";
  z_matcher_.DescribeNegationTo(out);
}

template<typename Scalar, typename Frame>
template<typename XMatcher, typename YMatcher, typename ZMatcher>
ComponentwiseMatcher3Impl<geometry::Bivector<Scalar, Frame> const&>::
ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                          YMatcher const& y_matcher,
                          ZMatcher const& z_matcher)
    : x_matcher_(SafeMatcherCast<Scalar>(x_matcher)),
      y_matcher_(SafeMatcherCast<Scalar>(y_matcher)),
      z_matcher_(SafeMatcherCast<Scalar>(z_matcher)) {}

template<typename Scalar, typename Frame>
bool ComponentwiseMatcher3Impl<geometry::Bivector<Scalar, Frame> const&>::
MatchAndExplain(geometry::Bivector<Scalar, Frame> const& actual,
                MatchResultListener* listener) const {
  bool const x_matches =
      x_matcher_.MatchAndExplain(actual.coordinates().x, listener);
  if (!x_matches) {
    *listener << " in the x coordinate; ";
  }
  bool const y_matches =
      y_matcher_.MatchAndExplain(actual.coordinates().y, listener);
  if (!y_matches) {
    *listener << " in the y coordinate; ";
  }
  bool const z_matches =
      z_matcher_.MatchAndExplain(actual.coordinates().z, listener);
  if (!z_matches) {
    *listener << " in the z coordinate; ";
  }
  return x_matches && y_matches && z_matches;
}

template<typename Scalar, typename Frame>
void ComponentwiseMatcher3Impl<geometry::Bivector<Scalar, Frame> const&>::
DescribeTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeTo(out);
  *out << " and y ";
  y_matcher_.DescribeTo(out);
  *out << " and z ";
  z_matcher_.DescribeTo(out);
}

template<typename Scalar, typename Frame>
void ComponentwiseMatcher3Impl<geometry::Bivector<Scalar, Frame> const&>::
DescribeNegationTo(std::ostream* out) const {
  *out << "x ";
  x_matcher_.DescribeNegationTo(out);
  *out << " or y ";
  y_matcher_.DescribeNegationTo(out);
  *out << " or z ";
  z_matcher_.DescribeNegationTo(out);
}

}  // namespace internal_componentwise
}  // namespace testing_utilities
}  // namespace principia
