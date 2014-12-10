#include "physics/degrees_of_freedom.hpp"

#include <vector>

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"

using principia::geometry::Displacement;
using principia::geometry::Position;
using principia::geometry::Velocity;
using principia::quantities::Entropy;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::SIUnit;
using principia::testing_utilities::AlmostEquals;

namespace principia {
namespace physics {

class DegreesOfFreedomTest : public testing::Test {
 protected:
  struct World;

  DegreesOfFreedomTest()
      : testing::Test(),
        d1_(origin_ + Displacement<World>({1 * SIUnit<Length>(),
                                           2 * SIUnit<Length>(),
                                           3 * SIUnit<Length>()}),
            Velocity<World>({10 * SIUnit<Speed>(),
                             20 * SIUnit<Speed>(),
                             30 * SIUnit<Speed>()})),
        d2_(origin_ + Displacement<World>({4 * SIUnit<Length>(),
                                           -5 * SIUnit<Length>(),
                                           6 * SIUnit<Length>()}),
            Velocity<World>({40 * SIUnit<Speed>(),
                             50 * SIUnit<Speed>(),
                             -60 * SIUnit<Speed>()})),
        d3_(origin_ + Displacement<World>({-7 * SIUnit<Length>(),
                                           8 * SIUnit<Length>(),
                                           -9 * SIUnit<Length>()}),
            Velocity<World>({-70 * SIUnit<Speed>(),
                             -80 * SIUnit<Speed>(),
                             -90 * SIUnit<Speed>()})) {}

  Position<World> origin_;
  DegreesOfFreedom<World> d1_;
  DegreesOfFreedom<World> d2_;
  DegreesOfFreedom<World> d3_;
};

TEST_F(DegreesOfFreedomTest, BarycentreError) {
  // The <> seem to confuse EXPECT_DEATH, hence the lambda.
  auto barycentre =
      [](std::vector<DegreesOfFreedom<World>> const& degrees_of_freedom,
         std::vector<Entropy> const& weights) -> DegreesOfFreedom<World> {
    return Barycentre<World, Entropy>(degrees_of_freedom, weights);
  };
  EXPECT_DEATH({
    barycentre({d1_, d2_, d3_}, {3 * SIUnit<Entropy>(), 4 * SIUnit<Entropy>()});
  }, "unequal sizes");
  EXPECT_DEATH({
    barycentre({}, {});
  }, "Empty input");
  EXPECT_DEATH({
    DegreesOfFreedom<World>::BarycentreCalculator<Entropy> calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
}

TEST_F(DegreesOfFreedomTest, Barycentre) {
  DegreesOfFreedom<World> const barycentre =
      Barycentre<World, Entropy>({d1_, d2_, d3_},
                                 {3 * SIUnit<Entropy>(),
                                  4 * SIUnit<Entropy>(),
                                  5 * SIUnit<Entropy>()});
  EXPECT_THAT(barycentre.position - origin_,
              AlmostEquals(
                  Displacement<World>({(-4.0 / 3.0) * SIUnit<Length>(),
                                       (13.0 / 6.0) * SIUnit<Length>(),
                                       -1.0 * SIUnit<Length>()}), 1));
  EXPECT_THAT(barycentre.velocity,
              AlmostEquals(
                  Velocity<World>({(-40.0 / 3.0) * SIUnit<Speed>(),
                                   (-35.0 / 3.0) * SIUnit<Speed>(),
                                   -50.0 * SIUnit<Speed>()}), 1));
}

TEST_F(DegreesOfFreedomTest, BarycentreCalculator) {
  DegreesOfFreedom<World>::BarycentreCalculator<double> calculator;
  calculator.Add(d1_, 3);
  DegreesOfFreedom<World> barycentre = calculator.Get();
  EXPECT_THAT(barycentre.position - origin_,
              AlmostEquals(d1_.position - origin_, 0));
  EXPECT_THAT(barycentre.velocity,
              AlmostEquals(d1_.velocity, 0));
  calculator.Add(d2_, 4);
  barycentre = calculator.Get();
  EXPECT_THAT(barycentre.position - origin_,
              AlmostEquals(
                  Displacement<World>({(19.0 / 7.0) * SIUnit<Length>(),
                                       -2.0 * SIUnit<Length>(),
                                       (33.0 / 7.0) * SIUnit<Length>()}), 0));
  EXPECT_THAT(barycentre.velocity,
              AlmostEquals(
                  Velocity<World>({(190.0 / 7.0) * SIUnit<Speed>(),
                                   (260.0 / 7.0) * SIUnit<Speed>(),
                                   (-150.0 / 7.0) * SIUnit<Speed>()}), 0));
  calculator.Add(d3_, 5);
  barycentre = calculator.Get();
  EXPECT_THAT(barycentre.position - origin_,
              AlmostEquals(
                  Displacement<World>({(-4.0 / 3.0) * SIUnit<Length>(),
                                        (13.0 / 6.0) * SIUnit<Length>(),
                                        -1.0 * SIUnit<Length>()}), 0));
  EXPECT_THAT(barycentre.velocity,
              AlmostEquals(
                  Velocity<World>({(-40.0 / 3.0) * SIUnit<Speed>(),
                                   (-35.0 / 3.0) * SIUnit<Speed>(),
                                   -50.0 * SIUnit<Speed>()}), 0));
}

}  // namespace physics
}  // namespace principia
