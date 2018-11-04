
#pragma once

#include <cfloat>
#include <cstdint>

#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/pair.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rp2_point.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_componentwise {

using geometry::R3Element;
using ::testing::Matcher;
using ::testing::MatcherInterface;
using ::testing::MatchResultListener;
using ::testing::PolymorphicMatcher;

template<typename T1Matcher, typename T2Matcher>
class ComponentwiseMatcher2;

template<typename XMatcher, typename YMatcher, typename ZMatcher>
class ComponentwiseMatcher3;

template<typename T1Matcher, typename T2Matcher>
PolymorphicMatcher<ComponentwiseMatcher2<T1Matcher, T2Matcher>>
Componentwise(T1Matcher const& t1_matcher,
              T2Matcher const& t2_matcher);

template<typename XMatcher, typename YMatcher, typename ZMatcher>
PolymorphicMatcher<ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher);

template<typename T1Matcher, typename T2Matcher>
class ComponentwiseMatcher2 final {
 public:
  ComponentwiseMatcher2(T1Matcher const& t1_matcher,
                        T2Matcher const& t2_matcher);

  template <typename PairType>
  operator Matcher<PairType> () const;

 private:
  T1Matcher const t1_matcher_;
  T2Matcher const t2_matcher_;
};

template<typename T1Matcher, typename T2Matcher, typename T3Matcher>
class ComponentwiseMatcher3 final {
 public:
  ComponentwiseMatcher3(T1Matcher const& t1_matcher,
                        T2Matcher const& t2_matcher,
                        T3Matcher const& t3_matcher);

  template <typename TripleType>
  operator Matcher<TripleType> () const;

 private:
  T1Matcher const t1_matcher_;
  T2Matcher const t2_matcher_;
  T3Matcher const t3_matcher_;
};

template<typename PairType>
class ComponentwiseMatcher2Impl;

template<typename T1, typename T2>
class ComponentwiseMatcher2Impl<geometry::Pair<T1, T2>> final
    : public MatcherInterface<geometry::Pair<T1, T2>> {
 public:
  template<typename T1Matcher, typename T2Matcher>
  ComponentwiseMatcher2Impl(T1Matcher const& t1_matcher,
                            T2Matcher const& t2_matcher);

  // Note that at this point this is only useful for vector/vector pairs as we
  // don't have matchers for |Point|.
  bool MatchAndExplain(geometry::Pair<T1, T2> const& actual,
                       MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  Matcher<T1> const t1_matcher_;
  Matcher<T2> const t2_matcher_;
};

template<typename Scalar, typename Frame>
class ComponentwiseMatcher2Impl<geometry::RP2Point<Scalar, Frame>> final
    : public MatcherInterface<geometry::RP2Point<Scalar, Frame>> {
 public:
  template<typename XMatcher, typename YMatcher>
  ComponentwiseMatcher2Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher);

  bool MatchAndExplain(geometry::RP2Point<Scalar, Frame> const& actual,
                       MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  Matcher<Scalar> const x_matcher_;
  Matcher<Scalar> const y_matcher_;
};

template<typename TripleType>
class ComponentwiseMatcher3Impl;

template<typename Scalar>
class ComponentwiseMatcher3Impl<geometry::R3Element<Scalar>> final
    : public MatcherInterface<geometry::R3Element<Scalar>> {
 public:
  template<typename XMatcher, typename YMatcher, typename ZMatcher>
  ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher,
                            ZMatcher const& z_matcher);

  bool MatchAndExplain(geometry::R3Element<Scalar> const& actual,
                       MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  Matcher<Scalar> const x_matcher_;
  Matcher<Scalar> const y_matcher_;
  Matcher<Scalar> const z_matcher_;
};

template<typename Scalar, typename Frame>
class ComponentwiseMatcher3Impl<geometry::Vector<Scalar, Frame>> final
    : public MatcherInterface<geometry::Vector<Scalar, Frame>> {
 public:
  template<typename XMatcher, typename YMatcher, typename ZMatcher>
  ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher,
                            ZMatcher const& z_matcher);

  bool MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                       MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  Matcher<Scalar> const x_matcher_;
  Matcher<Scalar> const y_matcher_;
  Matcher<Scalar> const z_matcher_;
};

template<typename Scalar, typename Frame>
class ComponentwiseMatcher3Impl<geometry::Bivector<Scalar, Frame>> final
    : public MatcherInterface<geometry::Bivector<Scalar, Frame>> {
 public:
  template<typename XMatcher, typename YMatcher, typename ZMatcher>
  ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher,
                            ZMatcher const& z_matcher);

  bool MatchAndExplain(geometry::Bivector<Scalar, Frame> const& actual,
                       MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  Matcher<Scalar> const x_matcher_;
  Matcher<Scalar> const y_matcher_;
  Matcher<Scalar> const z_matcher_;
};

}  // namespace internal_componentwise

using internal_componentwise::Componentwise;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/componentwise_body.hpp"
