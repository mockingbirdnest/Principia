
#include "physics/degrees_of_freedom.hpp"

#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
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
using geometry::Frame;
using geometry::Position;
using geometry::Velocity;
using quantities::Entropy;
using quantities::Length;
using quantities::Speed;
using testing_utilities::Componentwise;
using ::testing::Eq;
namespace si = quantities::si;

class DegreesOfFreedomTest : public testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  DegreesOfFreedomTest()
      : d1_(origin_ + Displacement<World>({1 * si::Unit<Length>,
                                           2 * si::Unit<Length>,
                                           3 * si::Unit<Length>}),
            Velocity<World>({10 * si::Unit<Speed>,
                             20 * si::Unit<Speed>,
                             30 * si::Unit<Speed>})),
        d2_(origin_ + Displacement<World>({4 * si::Unit<Length>,
                                           -5 * si::Unit<Length>,
                                           6 * si::Unit<Length>}),
            Velocity<World>({40 * si::Unit<Speed>,
                             50 * si::Unit<Speed>,
                             -60 * si::Unit<Speed>})),
        d3_(origin_ + Displacement<World>({-7 * si::Unit<Length>,
                                           8 * si::Unit<Length>,
                                           -9 * si::Unit<Length>}),
            Velocity<World>({-70 * si::Unit<Speed>,
                             -80 * si::Unit<Speed>,
                             -90 * si::Unit<Speed>})) {}

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
    barycentre({d1_, d2_, d3_}, {3 * si::Unit<Entropy>, 4 * si::Unit<Entropy>});
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
          {3 * si::Unit<Entropy>,
           4 * si::Unit<Entropy>,
           5 * si::Unit<Entropy>});
  EXPECT_THAT(barycentre,
              Componentwise(
                  Eq(origin_ +
                     Displacement<World>({(-4.0 / 3.0) * si::Unit<Length>,
                                          (13.0 / 6.0) * si::Unit<Length>,
                                          -1.0 * si::Unit<Length>})),
                  Eq(Velocity<World>({(-40.0 / 3.0) * si::Unit<Speed>,
                                      (-35.0 / 3.0) * si::Unit<Speed>,
                                      -50.0 * si::Unit<Speed>}))));
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
                     Displacement<World>({(19.0 / 7.0) * si::Unit<Length>,
                                          -2.0 * si::Unit<Length>,
                                          (33.0 / 7.0) * si::Unit<Length>})),
                  Eq(Velocity<World>({(190.0 / 7.0) * si::Unit<Speed>,
                                      (260.0 / 7.0) * si::Unit<Speed>,
                                      (-150.0 / 7.0) * si::Unit<Speed>}))));
  calculator.Add(d3_, 5);
  barycentre = calculator.Get();
  EXPECT_THAT(barycentre,
              Componentwise(
                  Eq(origin_ +
                     Displacement<World>({(-4.0 / 3.0) * si::Unit<Length>,
                                          (13.0 / 6.0) * si::Unit<Length>,
                                          -1.0 * si::Unit<Length>})),
                  Eq(Velocity<World>({(-40.0 / 3.0) * si::Unit<Speed>,
                                      (-35.0 / 3.0) * si::Unit<Speed>,
                                      -50.0 * si::Unit<Speed>}))));
}

}  // namespace internal_degrees_of_freedom
}  // namespace physics
}  // namespace principia
