#include "ksp_plugin/manœuvre.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {

using quantities::Pow;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;

namespace ksp_plugin {

class ManœuvreTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;
};

TEST_F(ManœuvreTest, UnitBurn) {
  Instant const t0 = Instant();
  Vector<double, World> e_y({0, 1, 0});
  Manœuvre<World> manœuvre(1 * Newton, 1 * Kilogram, 1 * Metre / Second, e_y);
  EXPECT_EQ(1 * Newton, manœuvre.thrust());
  EXPECT_EQ(1 * Kilogram, manœuvre.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre.effective_exhaust_velocity());
  EXPECT_EQ(e_y, manœuvre.direction());
  EXPECT_EQ(1 * Kilogram / Second, manœuvre.mass_flow());
  manœuvre.set_duration(1 * Second);
  EXPECT_EQ(1 * Second, manœuvre.duration());
  EXPECT_EQ(1 * Metre / Second, manœuvre.Δv());
  EXPECT_EQ(1 * Second, manœuvre.time_to_half_Δv());
  EXPECT_EQ(1 * Kilogram, manœuvre.final_mass());
  manœuvre.set_initial_time(t0);
  EXPECT_EQ(t0 + 1 * Second, manœuvre.final_time());
  EXPECT_EQ(t0 + 1 * Second, manœuvre.time_of_half_Δv());
  EXPECT_EQ(1 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.initial_time()).Norm());
  EXPECT_EQ(1 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.time_of_half_Δv()).Norm());
  EXPECT_EQ(1 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.final_time()).Norm());
}

}  // namespace ksp_plugin
}  // namespace principia
