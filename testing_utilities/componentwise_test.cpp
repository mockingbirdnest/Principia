
#include "testing_utilities/componentwise.hpp"

#include <limits>

#include "geometry/grassmann.hpp"
#include "geometry/pair.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rp2_point.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using geometry::Bivector;
using geometry::Pair;
using geometry::R3Element;
using geometry::RP2Point;
using geometry::Vector;
using quantities::Action;
using quantities::Amount;
using quantities::Length;
using quantities::Speed;
using quantities::SIUnit;
using quantities::si::Metre;
using quantities::si::Second;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Matcher;
using ::testing::Not;
using ::testing::_;

namespace testing_utilities {

struct World;

class ComponentwiseTest : public testing::Test {};

TEST_F(ComponentwiseTest, R3Element) {
  R3Element<double> r({1.0 + 1.0e-12, 1.0e-10, 3.5});
  EXPECT_THAT(r, Componentwise(AlmostEquals(1.0, 4504),
                               VanishesBefore(1.0, 450360),
                               Eq(3.5)));
  EXPECT_THAT(r, Componentwise(AlmostEquals(1.0, 4504),
                               VanishesBefore(1.0, 450360),
                               Not(Eq(2.5))));
  EXPECT_THAT(r, Not(Componentwise(AlmostEquals(1.0, 4),
                                   VanishesBefore(1.0, 4),
                                   Eq(2.5))));
}

TEST_F(ComponentwiseTest, Grassmann) {
  Vector<Length, World> v({(1.0 + 1.0e-12) * Metre,
                            1.0e-10 * Metre,
                            3.5 * Metre});
  EXPECT_THAT(v, Componentwise(AlmostEquals(1.0 * Metre, 4504),
                               VanishesBefore(1.0 * Metre, 450360),
                               Eq(3.5 * Metre)));
  EXPECT_THAT(v, Not(Componentwise(AlmostEquals(1.0 * Metre, 4),
                                   VanishesBefore(1.0 * Metre, 4),
                                   Eq(2.5 * Metre))));
  Bivector<Length, World> b({(1.0 + 1.0e-12) * Metre,
                              1.0e-10 * Metre,
                              3.5 * Metre});
  EXPECT_THAT(b, Componentwise(AlmostEquals(1.0 * Metre, 4504),
                               VanishesBefore(1.0 * Metre, 450360),
                               Eq(3.5 * Metre)));
  EXPECT_THAT(b, Not(Componentwise(AlmostEquals(1.0 * Metre, 4),
                                   VanishesBefore(1.0 * Metre, 4),
                                   Eq(2.5 * Metre))));
}

TEST_F(ComponentwiseTest, Pair) {
  using VV = Pair<Vector<Action, World>, Vector<Amount, World>>;
  VV vv(Vector<Action, World>({(1.0 + 1.0e-12) * SIUnit<Action>(),
                                1.0e-10 *  SIUnit<Action>(),
                                3.5 *  SIUnit<Action>()}),
        Vector<Amount, World>({(1.0 + 1.0e-12) * SIUnit<Amount>(),
                                (2.0 + 1.0e-10) *  SIUnit<Amount>(),
                                 3.5 *  SIUnit<Amount>()}));
  EXPECT_THAT(vv, Componentwise(
                      Componentwise(
                          AlmostEquals(1.0 * SIUnit<Action>(), 4504),
                          VanishesBefore(1.0 * SIUnit<Action>(), 450360),
                          Eq(3.5 * SIUnit<Action>())),
                      AlmostEquals(
                          Vector<Amount, World>({1.0 * SIUnit<Amount>(),
                                                 2.0 *  SIUnit<Amount>(),
                                                 3.5 *  SIUnit<Amount>()}),
                          225180)));
  EXPECT_THAT(vv, Not(Componentwise(
                      Componentwise(
                          AlmostEquals(1.0 * SIUnit<Action>(), 4504),
                          VanishesBefore(1.0 * SIUnit<Action>(), 450360),
                          Eq(2.5 * SIUnit<Action>())),
                      AlmostEquals(
                          Vector<Amount, World>({1.0 * SIUnit<Amount>(),
                                                 2.0 *  SIUnit<Amount>(),
                                                 3.5 *  SIUnit<Amount>()}),
                          2))));
}

TEST_F(ComponentwiseTest, Describe) {
  using RP2 = RP2Point<double, World>;
  using R3 = R3Element<double>;
  {
    std::ostringstream out;
    Matcher<R3 const&>(Componentwise(AlmostEquals(1.0, 2),
                              VanishesBefore(1.0, 4),
                              Eq(3.5))).DescribeTo(&out);
    EXPECT_EQ("x is within 2 to 2 ULPs of 1 and "
              "y vanishes before 1 to within 4 to 4 ULPs and "
              "z is equal to 3.5",
              out.str());
  }
  {
    std::ostringstream out;
    Matcher<R3 const&>(Componentwise(AlmostEquals(1.0, 2),
                              VanishesBefore(1.0, 4),
                              Eq(3.5))).DescribeNegationTo(&out);
    EXPECT_EQ("x is not within 2 to 2 ULPs of 1 or "
              "y does not vanish before 1 to within 4 to 4 ULP or "
              "z isn't equal to 3.5",
              out.str());
  }
  {
    std::ostringstream out;
    Matcher<RP2 const&>(Componentwise(AlmostEquals(1.0, 2),
                               VanishesBefore(1.0, 4))).DescribeTo(&out);
    EXPECT_EQ("t1 is within 2 to 2 ULPs of 1 and "
              "t2 vanishes before 1 to within 4 to 4 ULPs",
              out.str());
  }
  {
    std::ostringstream out;
    Matcher<RP2 const&>(Componentwise(AlmostEquals(1.0, 2),
                               VanishesBefore(1.0, 4)))
        .DescribeNegationTo(&out);
    EXPECT_EQ("t2 is not within 2 to 2 ULPs of 1 or "
              "t2 does not vanish before 1 to within 4 to 4 ULP",
              out.str());
  }
}

TEST_F(ComponentwiseTest, Variadic) {
  using VV = Pair<Vector<Length, World>, Vector<Speed, World>>;
  VV vv(Vector<Length, World>({1 * Metre,
                               2 * Metre,
                               3 * Metre}),
        Vector<Speed, World>({4 * Metre / Second,
                              5 * Metre / Second,
                              6 * Metre / Second}));
  EXPECT_THAT(
      vv,
      Componentwise(Componentwise(IsNear(1 * Metre),
                                  IsNear(2 * Metre),
                                  AllOf(Gt(2 * Metre), Lt(4 * Metre))),
                    Componentwise(IsNear(4 * Metre / Second),
                                  IsNear(5 * Metre / Second),
                                  AllOf(Gt(5 * Metre / Second),
                                        Lt(7 * Metre / Second)))));
}

TEST_F(ComponentwiseTest, Values) {
  using VV = Pair<Vector<Length, World>, Vector<Speed, World>>;
  VV vv(Vector<Length, World>({1 * Metre,
                               2 * Metre,
                               3 * Metre}),
        Vector<Speed, World>({4 * Metre / Second,
                              5 * Metre / Second,
                              6 * Metre / Second}));
  EXPECT_THAT(
      vv,
      Componentwise(Componentwise(1 * Metre,
                                  2 * Metre,
                                  3 * Metre),
                    Componentwise(4 * Metre / Second,
                                  5 * Metre / Second,
                                  6 * Metre / Second)));
}

TEST_F(ComponentwiseTest, Underscore) {
  using V = Vector<Length, World>;
  V v = Vector<Length, World>({1e-50 * Metre,
                               2e-50 * Metre,
                               3 * Metre});
  EXPECT_THAT(v,
              Componentwise(VanishesBefore(1 * Metre, 0),
                            VanishesBefore(1 * Metre, 0),
                            _));
}

}  // namespace testing_utilities
}  // namespace principia
