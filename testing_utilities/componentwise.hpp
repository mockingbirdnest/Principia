
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

template<typename T1Matcher, typename T2Matcher>
class ComponentwiseMatcher2;

template<typename XMatcher, typename YMatcher, typename ZMatcher>
class ComponentwiseMatcher3;

template<typename T1Matcher, typename T2Matcher>
testing::PolymorphicMatcher<ComponentwiseMatcher2<T1Matcher, T2Matcher>>
Componentwise(T1Matcher const& t1_matcher,
              T2Matcher const& t2_matcher);

template<typename XMatcher, typename YMatcher, typename ZMatcher>
testing::PolymorphicMatcher<ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher);

template<typename T1Matcher, typename T2Matcher>
class ComponentwiseMatcher2 final {
 public:
  explicit ComponentwiseMatcher2(T1Matcher const& t1_matcher,
                                 T2Matcher const& t2_matcher);

  // Note that at this point this is only useful for vector/vector pairs as we
  // don't have matchers for |Point|.
  template<typename T1, typename T2>
  bool MatchAndExplain(geometry::Pair<T1, T2> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::RP2Point<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  T1Matcher const& t1_matcher_;
  T2Matcher const& t2_matcher_;
};

template<typename XMatcher, typename YMatcher, typename ZMatcher>
class ComponentwiseMatcher3 final {
 public:
  explicit ComponentwiseMatcher3(XMatcher const& x_matcher,
                                 YMatcher const& y_matcher,
                                 ZMatcher const& z_matcher);

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

}  // namespace internal_componentwise

using internal_componentwise::Componentwise;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/componentwise_body.hpp"
