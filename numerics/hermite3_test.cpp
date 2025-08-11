#include "numerics/hermite3.hpp"

#include <utility>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace numerics {

using ::testing::ElementsAre;
using ::testing::Eq;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;

class Hermite3Test : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag, Inertial>;

  Instant const t0_;
};

TEST_F(Hermite3Test, Precomputed) {
  Hermite3<Length, Instant> h({t0_ + 1 * Second, t0_ + 2 * Second},
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

TEST_F(Hermite3Test, Quadratic) {
  Hermite3<Length, Instant> near_quadratic(
      {t0_ + 1 * Second, t0_ + 2 * Second},
      {0x1p-53 * Metre, 0 * Metre},
      {-1 * Metre / Second, 1 * Metre / Second});
  // These are the correctly-rounded extrema.
  EXPECT_THAT(
      near_quadratic.FindExtrema(),
      ElementsAre(t0_ - 0x1.5555555555552p51 * Second, t0_ + 1.5 * Second));
  Hermite3<Length, Instant> quadratic(
      {t0_ + 1 * Second, t0_ + 2 * Second},
      {0 * Metre, 0 * Metre},
      {-1 * Metre / Second, 1 * Metre / Second});
  EXPECT_THAT(quadratic.FindExtrema(),
              ElementsAre(InfinitePast, t0_ + 1.5 * Second));
}

TEST_F(Hermite3Test, Typed) {
  // Just here to check that the types work in the presence of affine spaces.
  Hermite3<Position<World>, Instant> h({t0_ + 1 * Second, t0_ + 2 * Second},
                                       {World::origin, World::origin},
                                       {World::unmoving, World::unmoving});

  EXPECT_EQ(World::origin, h.Evaluate(t0_ + 1.3 * Second));
  EXPECT_EQ(Velocity<World>(), h.EvaluateDerivative(t0_ + 1.7 * Second));
}

// A test that shows a poorly-conditioned case where the computed position at
// the end of the interpolation interval is different from the expected one.
// In particular, the computed position is negative at both ends.
TEST_F(Hermite3Test, Conditioning) {
  Hermite3<Length, Instant> h(
      {t0_ + 19418861.806896236 * Second, t0_ + 19418869.842261545 * Second},
      {-2.1383610158805017e-12 * Metre, 0 * Metre},
      {2.3308208035605881e-12 * Metre / Second,
       -1.6875631840376598e-12 * Metre / Second});

  EXPECT_EQ(-2.1383610158805017e-12 * Metre,
            h.Evaluate(t0_ + 19418861.806896236 * Second));
  EXPECT_GT(0 * Metre, h.Evaluate(t0_ + 19418869.842261545 * Second));
  EXPECT_EQ(2.3308208035605881e-12 * Metre / Second,
            h.EvaluateDerivative(t0_ + 19418861.806896236 * Second));
}

TEST_F(Hermite3Test, OneDimensionalInterpolationError) {
  std::vector<std::pair<double, double>> samples;
  for (double i = 2; i < 10; i += 1) {
    samples.push_back({1 / i, Pow<4>(1 / i)});
  }
  const auto not_a_quartic =
      Hermite3<double, double>(/*arguments=*/{0, 1},
                               /*values=*/{0, 1},
                               /*derivatives=*/{0, 4});
  // `not_a_quartic` has a root at 1/2, where the error is maximal.
  EXPECT_THAT(not_a_quartic.LInfinityError(
      samples,
      /*get_argument=*/[](auto&& pair) -> auto&& { return pair.first; },
      /*get_value=*/[](auto&& pair) -> auto&& { return pair.second; }),
      Eq(1 / 16.0));

  EXPECT_TRUE(not_a_quartic.LInfinityErrorIsWithin(
      samples,
      /*get_argument=*/[](auto&& pair) -> auto&& { return pair.first; },
      /*get_value=*/[](auto&& pair) -> auto&& { return pair.second; },
      /*tolerance=*/0.1));
  EXPECT_FALSE(not_a_quartic.LInfinityErrorIsWithin(
      samples,
      /*get_argument=*/[](auto&& pair) -> auto&& { return pair.first; },
      /*get_value=*/[](auto&& pair) -> auto&& { return pair.second; },
      /*tolerance=*/0.05));
}

TEST_F(Hermite3Test, ThreeDimensionalInterpolationError) {
  std::vector<std::pair<Instant, Position<World>>> samples;
  Instant const tmax = t0_ + π / 2 * Second;
  AngularFrequency const ω = 1 * Radian / Second;
  for (Instant t = t0_; t <= tmax; t += 1 / 32.0 * Second) {
    samples.push_back(
        {t, World::origin + Displacement<World>({Cos(ω * (t - t0_)) * Metre,
                                                 Sin(ω * (t - t0_)) * Metre,
                                                 0 * Metre})});
  }
  const auto not_a_circle = Hermite3<Position<World>, Instant>(
      /*arguments=*/{t0_, tmax},
      /*values=*/
      {World::origin + Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}),
       World::origin + Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre})},
      /*derivatives=*/
      {Velocity<World>(
           {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second}),
       Velocity<World>(
           {-1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second})});
  EXPECT_THAT(
      not_a_circle.LInfinityError(
          samples,
          /*get_argument=*/[](auto&& pair) -> auto&& { return pair.first; },
          /*get_value=*/[](auto&& pair) -> auto&& { return pair.second; }),
      IsNear(1.5_(1) * Centi(Metre)));

  EXPECT_TRUE(not_a_circle.LInfinityErrorIsWithin(
      samples,
      /*get_argument=*/[](auto&& pair) -> auto&& { return pair.first; },
      /*get_value=*/[](auto&& pair) -> auto&& { return pair.second; },
      /*tolerance=*/2 * Centi(Metre)));
  EXPECT_FALSE(not_a_circle.LInfinityErrorIsWithin(
      samples,
      /*get_argument=*/[](auto&& pair) -> auto&& { return pair.first; },
      /*get_value=*/[](auto&& pair) -> auto&& { return pair.second; },
      /*tolerance=*/1 * Centi(Metre)));
}

}  // namespace numerics
}  // namespace principia
