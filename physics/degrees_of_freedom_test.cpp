
#include "physics/degrees_of_freedom.hpp"

#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace physics {
namespace internal_degrees_of_freedom {

using geometry::Barycentre;
using geometry::BarycentreCalculator;
using geometry::Displacement;
using geometry::Position;
using geometry::Velocity;
using quantities::Entropy;
using quantities::Length;
using quantities::Speed;
using quantities::SIUnit;
using testing_utilities::Componentwise;
using ::testing::Eq;

class DegreesOfFreedomTest : public testing::Test {
 protected:
  struct World;

  DegreesOfFreedomTest()
      : d1_(origin_ + Displacement<World>({1 * SIUnit<Length>(),
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

using DegreesOfFreedomDeathTest = DegreesOfFreedomTest;

TEST_F(DegreesOfFreedomDeathTest, BarycentreError) {
  // The <> seem to confuse EXPECT_DEATH, hence the lambda.
  auto barycentre =
      [](std::vector<DegreesOfFreedom<World>> const& degrees_of_freedom,
         std::vector<Entropy> const& weights) -> DegreesOfFreedom<World> {
    return Barycentre<DegreesOfFreedom<World>, Entropy, std::vector>(
               degrees_of_freedom, weights);
  };
  EXPECT_DEATH({
    barycentre({d1_, d2_, d3_}, {3 * SIUnit<Entropy>(), 4 * SIUnit<Entropy>()});
  }, "unequal sizes");
  EXPECT_DEATH({
    barycentre({}, {});
  }, "Empty input");
  using DegreesOfFreedomBarycentreCalculator =
      BarycentreCalculator<DegreesOfFreedom<World>, Entropy>;
  EXPECT_DEATH({
    DegreesOfFreedomBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
}

TEST_F(DegreesOfFreedomTest, Output) {\
  EXPECT_EQ(DebugString(d1_),
    "{{+1.00000000000000000e+00 m, "
    "+2.00000000000000000e+00 m, "
    "+3.00000000000000000e+00 m}, "
    "{+1.00000000000000000e+01 m s^-1, "
    "+2.00000000000000000e+01 m s^-1, "
    "+3.00000000000000000e+01 m s^-1}}");
  RelativeDegreesOfFreedom<World> relative_degrees_of_freedom = d1_ - d2_;
  EXPECT_EQ(DebugString(relative_degrees_of_freedom),
    "{{-3.00000000000000000e+00 m, "
    "+7.00000000000000000e+00 m, "
    "-3.00000000000000000e+00 m}, "
    "{-3.00000000000000000e+01 m s^-1, "
    "-3.00000000000000000e+01 m s^-1, "
    "+9.00000000000000000e+01 m s^-1}}");
}

TEST_F(DegreesOfFreedomTest, Barycentre) {
  DegreesOfFreedom<World> const barycentre =
      Barycentre<DegreesOfFreedom<World>, Entropy, std::vector>(
          {d1_, d2_, d3_},
          {3 * SIUnit<Entropy>(),
           4 * SIUnit<Entropy>(),
           5 * SIUnit<Entropy>()});
  EXPECT_THAT(barycentre,
              Componentwise(
                  Eq(origin_ +
                     Displacement<World>({(-4.0 / 3.0) * SIUnit<Length>(),
                                          (13.0 / 6.0) * SIUnit<Length>(),
                                          -1.0 * SIUnit<Length>()})),
                  Eq(Velocity<World>({(-40.0 / 3.0) * SIUnit<Speed>(),
                                      (-35.0 / 3.0) * SIUnit<Speed>(),
                                      -50.0 * SIUnit<Speed>()}))));
}

TEST_F(DegreesOfFreedomTest, BarycentreCalculator) {
  BarycentreCalculator<DegreesOfFreedom<World>, double> calculator;
  calculator.Add(d1_, 3);
  DegreesOfFreedom<World> barycentre = calculator.Get();
  EXPECT_THAT(barycentre, Eq(d1_));
  calculator.Add(d2_, 4);
  barycentre = calculator.Get();
  EXPECT_THAT(barycentre,
              Componentwise(
                  Eq(origin_ +
                     Displacement<World>({(19.0 / 7.0) * SIUnit<Length>(),
                                          -2.0 * SIUnit<Length>(),
                                          (33.0 / 7.0) * SIUnit<Length>()})),
                  Eq(Velocity<World>({(190.0 / 7.0) * SIUnit<Speed>(),
                                      (260.0 / 7.0) * SIUnit<Speed>(),
                                      (-150.0 / 7.0) * SIUnit<Speed>()}))));
  calculator.Add(d3_, 5);
  barycentre = calculator.Get();
  EXPECT_THAT(barycentre,
              Componentwise(
                  Eq(origin_ +
                     Displacement<World>({(-4.0 / 3.0) * SIUnit<Length>(),
                                          (13.0 / 6.0) * SIUnit<Length>(),
                                          -1.0 * SIUnit<Length>()})),
                  Eq(Velocity<World>({(-40.0 / 3.0) * SIUnit<Speed>(),
                                      (-35.0 / 3.0) * SIUnit<Speed>(),
                                      -50.0 * SIUnit<Speed>()}))));
}

}  // namespace internal_degrees_of_freedom
}  // namespace physics
}  // namespace principia
