
#include "numerics/hermite3.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using geometry::Frame;
using geometry::Instant;
using geometry::Displacement;
using geometry::Position;
using geometry::Velocity;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::ElementsAre;

namespace numerics {

class Hermite3Test : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  Instant const t0_;
};

TEST_F(Hermite3Test, Precomputed) {
  Hermite3<Instant, Length> h({t0_ + 1 * Second, t0_ + 2 * Second},
                              {33 * Metre, 40 * Metre},
                              {-5 * Metre / Second, 6 * Metre / Second});

  EXPECT_EQ(33 * Metre, h.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(33.109375 * Metre, h.Evaluate(t0_ + 1.25 * Second));
  EXPECT_EQ(35.125 * Metre, h.Evaluate(t0_ + 1.5 * Second));
  EXPECT_EQ(37.828125 * Metre, h.Evaluate(t0_ + 1.75 * Second));
  EXPECT_EQ(40 * Metre, h.Evaluate(t0_ + 2 * Second));

  EXPECT_EQ(-5 * Metre / Second, h.EvaluateDerivative(t0_ + 1 * Second));
  EXPECT_EQ(5.0625 * Metre / Second, h.EvaluateDerivative(t0_ + 1.25 * Second));
  EXPECT_EQ(10.25 * Metre / Second, h.EvaluateDerivative(t0_ + 1.5 * Second));
  EXPECT_EQ(10.5625 * Metre / Second,
            h.EvaluateDerivative(t0_ + 1.75 * Second));
  EXPECT_EQ(6 * Metre / Second, h.EvaluateDerivative(t0_ + 2 * Second));

  EXPECT_THAT(h.FindExtrema(),
              ElementsAre(t0_ + ((64.0 - sqrt(430.0)) / 39.0) * Second,
                          t0_ + ((64.0 + sqrt(430.0)) / 39.0) * Second));
}

TEST_F(Hermite3Test, Typed) {
  // Just here to check that the types work in the presence of affine spaces.
  Hermite3<Instant, Position<World>> h({t0_ + 1 * Second, t0_ + 2 * Second},
                                       {World::origin, World::origin},
                                       {Velocity<World>(), Velocity<World>()});

  EXPECT_EQ(World::origin, h.Evaluate(t0_ + 1.3 * Second));
  EXPECT_EQ(Velocity<World>(), h.EvaluateDerivative(t0_ + 1.7 * Second));
}

// A test that shows a poorly-conditioned case where the computed position at
// the end of the interpolation interval is different from the expected one.
TEST_F(Hermite3Test, Conditioning) {
  Hermite3<Instant, Position<World>> h(
      {t0_ + 19418861.806896236 * Second, t0_ + 19418869.842261545 * Second},
      {World::origin + Displacement<World>({3.1953689924786266e-12 * Metre,
                                            6.1771019395071032e-12 * Metre,
                                            -2.1383610158805017e-12 * Metre}),
       World::origin},
      {Velocity<World>({1.2922104075254831e-12 * Metre / Second,
                        1.7453515790397885e-12 * Metre / Second,
                        2.3308208035605881e-12 * Metre / Second}),
       Velocity<World>({-4.0937095947499385e-12 * Metre / Second,
                        1.3971639017653770e-12 * Metre / Second,
                        -1.6875631840376598e-12 * Metre / Second})});

  EXPECT_EQ(
      World::origin + Displacement<World>({3.1953689924786266e-12 * Metre,
                                           6.1771019395071032e-12 * Metre,
                                           -2.1383610158805017e-12 * Metre}),
      h.Evaluate(t0_ + 19418861.806896236 * Second));
  EXPECT_NE(World::origin, h.Evaluate(t0_ + 19418869.842261545 * Second));
  EXPECT_EQ(Velocity<World>({1.2922104075254831e-12 * Metre / Second,
                             1.7453515790397885e-12 * Metre / Second,
                             2.3308208035605881e-12 * Metre / Second}),
            h.EvaluateDerivative(t0_ + 19418861.806896236 * Second));
  EXPECT_NE(Velocity<World>({-4.0937095947499385e-12 * Metre / Second,
                             1.3971639017653770e-12 * Metre / Second,
                             -1.6875631840376598e-12 * Metre / Second}),
            h.EvaluateDerivative(t0_ + 19418869.842261545 * Second));
}

}  // namespace numerics
}  // namespace principia
