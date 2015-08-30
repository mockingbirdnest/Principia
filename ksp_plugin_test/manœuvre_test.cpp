#include "ksp_plugin/manœuvre.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using quantities::Pow;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;
using testing_utilities::AlmostEquals;

namespace ksp_plugin {

class ManœuvreTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;
};

TEST_F(ManœuvreTest, TimedBurn) {
  Instant const t0 = Instant();
  Vector<double, World> e_y({0, 1, 0});
  Manœuvre<World> manœuvre(1 * Newton, 2 * Kilogram, 1 * Metre / Second, e_y);
  EXPECT_EQ(1 * Newton, manœuvre.thrust());
  EXPECT_EQ(2 * Kilogram, manœuvre.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre.effective_exhaust_velocity());
  EXPECT_EQ(e_y, manœuvre.direction());
  EXPECT_EQ(1 * Kilogram / Second, manœuvre.mass_flow());
  manœuvre.set_duration(1 * Second);
  EXPECT_EQ(1 * Second, manœuvre.duration());
  EXPECT_EQ(std::log(2) * Metre / Second, manœuvre.Δv());
  EXPECT_EQ((2 - Sqrt(2)) * Second, manœuvre.time_to_half_Δv());
  EXPECT_EQ(1 * Kilogram, manœuvre.final_mass());
  manœuvre.set_initial_time(t0);
  EXPECT_EQ(t0 + 1 * Second, manœuvre.final_time());
  EXPECT_EQ(t0 + (2 - Sqrt(2)) * Second, manœuvre.time_of_half_Δv());
  EXPECT_EQ(
      0 * Metre / Pow<2>(Second),
      manœuvre.acceleration()(manœuvre.initial_time() - 1 * Second).Norm());
  EXPECT_EQ(0.5 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.initial_time()).Norm());
  EXPECT_THAT(manœuvre.acceleration()(manœuvre.time_of_half_Δv()).Norm(),
              AlmostEquals(Sqrt(0.5) * Metre / Pow<2>(Second), 1));
  EXPECT_EQ(1 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.final_time()).Norm());
  EXPECT_EQ(0 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.final_time() + 1 * Second).Norm());
}

TEST_F(ManœuvreTest, TargetΔv) {
  Instant const t0 = Instant();
  Vector<double, World> e_y({0, 1, 0});
  Manœuvre<World> manœuvre(1 * Newton, 2 * Kilogram, 1 * Metre / Second, e_y);
  EXPECT_EQ(1 * Newton, manœuvre.thrust());
  EXPECT_EQ(2 * Kilogram, manœuvre.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre.effective_exhaust_velocity());
  EXPECT_EQ(e_y, manœuvre.direction());
  EXPECT_EQ(1 * Kilogram / Second, manœuvre.mass_flow());
  manœuvre.set_Δv(1 * Metre / Second);
  EXPECT_EQ(1 * Metre / Second, manœuvre.Δv());
  EXPECT_EQ((2 - 2 / e) * Second, manœuvre.duration());
  EXPECT_EQ((2 - 2 / Sqrt(e)) * Second, manœuvre.time_to_half_Δv());
  EXPECT_EQ((2 / e) * Kilogram, manœuvre.final_mass());
  manœuvre.set_initial_time(t0);
  EXPECT_EQ(t0 + (2 - 2 / e) * Second, manœuvre.final_time());
  EXPECT_EQ(t0 + (2 - 2 / Sqrt(e)) * Second, manœuvre.time_of_half_Δv());
  EXPECT_EQ(
      0 * Metre / Pow<2>(Second),
      manœuvre.acceleration()(manœuvre.initial_time() - 1 * Second).Norm());
  EXPECT_EQ(0.5 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.initial_time()).Norm());
  EXPECT_THAT(manœuvre.acceleration()(manœuvre.time_of_half_Δv()).Norm(),
              AlmostEquals(Sqrt(0.5) * Metre / Pow<2>(Second), 1));
  EXPECT_EQ((e / 2) * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.final_time()).Norm());
  EXPECT_EQ(0 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.final_time() + 1 * Second).Norm());
}

}  // namespace ksp_plugin
}  // namespace principia
