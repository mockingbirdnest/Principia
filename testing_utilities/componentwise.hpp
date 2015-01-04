#pragma once

#include <cfloat>
#include <cstdint>

#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/pair.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"

using principia::geometry::R3Element;

namespace principia {
namespace testing_utilities {

template<typename QMatcher, typename PMatcher>
class ComponentwiseMatcher2;

template<typename XMatcher, typename YMatcher, typename ZMatcher>
class ComponentwiseMatcher3;

template<typename QMatcher, typename PMatcher>
testing::PolymorphicMatcher<ComponentwiseMatcher2<QMatcher, PMatcher>>
Componentwise(QMatcher const& q_matcher,
              PMatcher const& p_matcher);

template<typename XMatcher, typename YMatcher, typename ZMatcher>
testing::PolymorphicMatcher<ComponentwiseMatcher3<XMatcher, YMatcher, ZMatcher>>
Componentwise(XMatcher const& x_matcher,
              YMatcher const& y_matcher,
              ZMatcher const& z_matcher);

template<typename QMatcher, typename PMatcher>
class ComponentwiseMatcher2 {
 public:
  explicit ComponentwiseMatcher2(QMatcher const& q_matcher,
                                 PMatcher const& p_matcher);
  ~ComponentwiseMatcher2() = default;

  template<typename T1, typename T2>
  bool MatchAndExplain(geometry::Pair<T1, T2> const& actual,
                       testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  QMatcher const& q_matcher_;
  PMatcher const& p_matcher_;
};

template<typename XMatcher, typename YMatcher, typename ZMatcher>
class ComponentwiseMatcher3 {
 public:
  explicit ComponentwiseMatcher3(XMatcher const& x_matcher,
                                 YMatcher const& y_matcher,
                                 ZMatcher const& z_matcher);
  ~ComponentwiseMatcher3() = default;

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

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/componentwise_body.hpp"
