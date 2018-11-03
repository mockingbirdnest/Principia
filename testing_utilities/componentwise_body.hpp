
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
using ::testing::Matcher;

// In order to call |Describe...To| on the various matchers we need to convert
// them to some |Matcher<T>|.  However, it is quite difficult to figure out
// what |T| is because the class |ComponentwiseMatcher| and the factory
// |Componentwise| cannot take it as a parameter (the template deduction would
// fail and the usages would get very ugly).  To make things worse, the actual
// type of the matcher depends on whether it is polymorphic, monomorphic or
// some other internal helper.  We obtain |T| by peeling away the layers of
// templates around it.
// TODO(phl): This is horribly complicated.  It is also very brittle as it
// depends on the exact structure of templates used in gmock-matchers.h.  One
// day I'll understand how the matchers work.

template<typename T>
struct MatcherParameterType : not_constructible {
  using type = T;
};

template<typename T, template<typename> class U>
struct MatcherParameterType<U<T>> : not_constructible {
  using type = typename MatcherParameterType<T>::type;
};

template<typename T1, typename T2, template<typename, typename> class U>
struct MatcherParameterType<U<T1, T2>> : not_constructible {
  // TODO(phl): Why T1?
  using type = typename MatcherParameterType<T1>::type;
};

template<
  template<template<typename...> class, typename, typename...> class U,
  template<typename> class V, typename T1, typename... T2>
struct MatcherParameterType<U<V, T1, T2...>> : not_constructible {
  // TODO(phl): Why T1?
  using type = typename MatcherParameterType<V<T1>>::type;
};

// |type| must be a type for which we implement MatchAndExplain.  We don't care
// which one exactly, since |MatcherParameterType| is only used for describing
// the matchers.  So we pick the simplest, |R3Element|.
template<typename XMatcher, typename YMatcher, typename ZMatcher>
struct MatcherParameterType<ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>>
    : not_constructible {
  using type =
      geometry::R3Element<typename MatcherParameterType<XMatcher>::type>;
};

// And now the cases that we *don't* want to peel away.  Yes, this smells a bit.
template<typename T>
struct MatcherParameterType<Point<T>> : not_constructible {
  using type = Point<T>;
};

template<typename T>
struct MatcherParameterType<Quantity<T>> : not_constructible {
  using type = Quantity<T>;
};

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
ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>::ComponentwiseMatcher3(
    XMatcher const& x_matcher,
    YMatcher const& y_matcher,
    ZMatcher const& z_matcher)
    : x_matcher_(x_matcher),
      y_matcher_(y_matcher),
      z_matcher_(z_matcher) {}

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
