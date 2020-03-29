
#pragma once

#include <cfloat>
#include <cstdint>

#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/pair.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rp2_point.hpp"
#include "gmock/gmock.h"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_componentwise {

using geometry::Displacement;
using geometry::Position;
using geometry::R3Element;
using geometry::Velocity;
using ::testing::Matcher;
using ::testing::MatcherInterface;
using ::testing::MatchResultListener;

template<typename T1Matcher, typename T2Matcher>
class ComponentwiseMatcher2;

template<typename T1Matcher, typename T2Matcher, typename T3Matcher>
class ComponentwiseMatcher3;

template<typename T1Matcher, typename T2Matcher>
ComponentwiseMatcher2<T1Matcher, T2Matcher>
Componentwise(T1Matcher const& t1_matcher,
              T2Matcher const& t2_matcher);

template<typename XMatcher, typename YMatcher, typename ZMatcher>
ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher);

template<typename T1Matcher, typename T2Matcher>
class ComponentwiseMatcher2 final {
 public:
  ComponentwiseMatcher2(T1Matcher t1_matcher,
                        T2Matcher t2_matcher);

  template <typename PairType>
  operator Matcher<PairType> () const;

 private:
  T1Matcher const t1_matcher_;
  T2Matcher const t2_matcher_;
};

template<typename T1Matcher, typename T2Matcher, typename T3Matcher>
class ComponentwiseMatcher3 final {
 public:
  ComponentwiseMatcher3(T1Matcher t1_matcher,
                        T2Matcher t2_matcher,
                        T3Matcher t3_matcher);

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
class ComponentwiseMatcher2Impl<geometry::Pair<T1, T2> const&> final
    : public MatcherInterface<geometry::Pair<T1, T2> const&> {
 public:
  template<typename T1Matcher, typename T2Matcher>
  ComponentwiseMatcher2Impl(T1Matcher const& t1_matcher,
                            T2Matcher const& t2_matcher);

  bool MatchAndExplain(geometry::Pair<T1, T2> const& actual,
                       MatchResultListener* listener) const override;

  void DescribeTo(std::ostream* out) const override;
  void DescribeNegationTo(std::ostream* out) const override;

 private:
  Matcher<T1> const t1_matcher_;
  Matcher<T2> const t2_matcher_;
};

template<typename Frame>
class ComponentwiseMatcher2Impl<physics::DegreesOfFreedom<Frame> const&> final
    : public MatcherInterface<physics::DegreesOfFreedom<Frame> const&> {
 public:
  template<typename QMatcher, typename PMatcher>
  ComponentwiseMatcher2Impl(QMatcher const& q_matcher,
                            PMatcher const& p_matcher);

  bool MatchAndExplain(physics::DegreesOfFreedom<Frame> const& actual,
                       MatchResultListener* listener) const override;

  void DescribeTo(std::ostream* out) const override;
  void DescribeNegationTo(std::ostream* out) const override;

 private:
  Matcher<Position<Frame>> const q_matcher_;
  Matcher<Velocity<Frame>> const p_matcher_;
};

template<typename Frame>
class ComponentwiseMatcher2Impl<physics::RelativeDegreesOfFreedom<Frame> const&>
    final
    : public MatcherInterface<physics::RelativeDegreesOfFreedom<Frame> const&> {
 public:
  template<typename QMatcher, typename PMatcher>
  ComponentwiseMatcher2Impl(QMatcher const& q_matcher,
                            PMatcher const& p_matcher);

  bool MatchAndExplain(physics::RelativeDegreesOfFreedom<Frame> const& actual,
                       MatchResultListener* listener) const override;

  void DescribeTo(std::ostream* out) const override;
  void DescribeNegationTo(std::ostream* out) const override;

 private:
  Matcher<Displacement<Frame>> const q_matcher_;
  Matcher<Velocity<Frame>> const p_matcher_;
};

template<typename Scalar, typename Frame>
class ComponentwiseMatcher2Impl<geometry::RP2Point<Scalar, Frame> const&> final
    : public MatcherInterface<geometry::RP2Point<Scalar, Frame> const&> {
 public:
  template<typename XMatcher, typename YMatcher>
  ComponentwiseMatcher2Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher);

  bool MatchAndExplain(geometry::RP2Point<Scalar, Frame> const& actual,
                       MatchResultListener* listener) const override;

  void DescribeTo(std::ostream* out) const override;
  void DescribeNegationTo(std::ostream* out) const override;

 private:
  Matcher<Scalar> const x_matcher_;
  Matcher<Scalar> const y_matcher_;
};

template<typename TripleType>
class ComponentwiseMatcher3Impl;

template<typename Scalar>
class ComponentwiseMatcher3Impl<geometry::R3Element<Scalar> const&> final
    : public MatcherInterface<geometry::R3Element<Scalar> const&> {
 public:
  template<typename XMatcher, typename YMatcher, typename ZMatcher>
  ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher,
                            ZMatcher const& z_matcher);

  bool MatchAndExplain(geometry::R3Element<Scalar> const& actual,
                       MatchResultListener* listener) const override;

  void DescribeTo(std::ostream* out) const override;
  void DescribeNegationTo(std::ostream* out) const override;

 private:
  Matcher<Scalar> const x_matcher_;
  Matcher<Scalar> const y_matcher_;
  Matcher<Scalar> const z_matcher_;
};

template<typename Scalar, typename Frame>
class ComponentwiseMatcher3Impl<geometry::Vector<Scalar, Frame> const&> final
    : public MatcherInterface<geometry::Vector<Scalar, Frame> const&> {
 public:
  template<typename XMatcher, typename YMatcher, typename ZMatcher>
  ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher,
                            ZMatcher const& z_matcher);

  bool MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                       MatchResultListener* listener) const override;

  void DescribeTo(std::ostream* out) const override;
  void DescribeNegationTo(std::ostream* out) const override;

 private:
  Matcher<Scalar> const x_matcher_;
  Matcher<Scalar> const y_matcher_;
  Matcher<Scalar> const z_matcher_;
};

template<typename Scalar, typename Frame>
class ComponentwiseMatcher3Impl<geometry::Bivector<Scalar, Frame> const&> final
    : public MatcherInterface<geometry::Bivector<Scalar, Frame> const&> {
 public:
  template<typename XMatcher, typename YMatcher, typename ZMatcher>
  ComponentwiseMatcher3Impl(XMatcher const& x_matcher,
                            YMatcher const& y_matcher,
                            ZMatcher const& z_matcher);

  bool MatchAndExplain(geometry::Bivector<Scalar, Frame> const& actual,
                       MatchResultListener* listener) const override;

  void DescribeTo(std::ostream* out) const override;
  void DescribeNegationTo(std::ostream* out) const override;

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
