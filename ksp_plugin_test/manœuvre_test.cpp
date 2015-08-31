#include "ksp_plugin/manœuvre.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using quantities::Pow;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;
using uk::Pound;
using uk::PoundForce;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

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
  EXPECT_EQ(t0, manœuvre.initial_time());
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

  manœuvre.set_time_of_half_Δv(t0);
  EXPECT_EQ(t0 - (2 - 2 / Sqrt(e)) * Second, manœuvre.initial_time());
  EXPECT_EQ(t0 + (2 / Sqrt(e) - 2 / e) * Second, manœuvre.final_time());
  EXPECT_EQ(t0, manœuvre.time_of_half_Δv());
  EXPECT_EQ(
      0 * Metre / Pow<2>(Second),
      manœuvre.acceleration()(manœuvre.initial_time() - 1 * Second).Norm());
  EXPECT_EQ(0.5 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.initial_time()).Norm());
  EXPECT_EQ(Sqrt(e) / 2 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.time_of_half_Δv()).Norm());
  EXPECT_EQ((e / 2) * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.final_time()).Norm());
  EXPECT_EQ(0 * Metre / Pow<2>(Second),
            manœuvre.acceleration()(manœuvre.final_time() + 1 * Second).Norm());
}

TEST_F(ManœuvreTest, Apollo8SIVB) {
  // Data from NASA's Saturn V Launch Vehicle, Flight Evaluation Report AS-503,
  // Apollo 8 Mission (1969),
  // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690015314.pdf.
  // We use the reconstructed or actual values.

  // Table 2-2. Significant Event Times Summary.
  Instant const range_zero;
  Instant const s_ivb_1st_ignition          = range_zero +    528.29 * Second;
  Instant const s_ivb_1st_90_percent_thrust = range_zero +    530.53 * Second;
  Instant const s_ivb_1st_eco               = range_zero +    684.98 * Second;
  // Inititiate S-IVB Restart Sequence and Start of Time Base 6 (T6).
  Instant const t6                          = range_zero +   9659.54 * Second;
  Instant const s_ivb_2nd_ignition          = range_zero + 10'237.79 * Second;
  Instant const s_ivb_2nd_90_percent_thrust = range_zero + 10'240.02 * Second;
  Instant const s_ivb_2nd_eco               = range_zero + 10'555.51 * Second;

  // From Table 7-2. S-IVB Steady State Performance - First Burn.
  Force thrust_1st                  = 901'557 * Newton;
  Speed specific_impulse_1st        = 4'204.1 * Newton * Second / Kilogram;
  Variation<Mass> lox_flowrate_1st  = 178.16 * Kilogram / Second;
  Variation<Mass> fuel_flowrate_1st = 36.30 * Kilogram / Second;

  // From Table 7-7. S-IVB Steady State Performance - Second Burn.
  Force thrust_2nd                  = 897'548 * Newton;
  Speed specific_impulse_2nd        = 4199.2 * Newton * Second / Kilogram;
  Variation<Mass> lox_flowrate_2nd  = 177.70 * Kilogram / Second;
  Variation<Mass> fuel_flowrate_2nd = 36.01 * Kilogram / Second;

  // Table 21-5. Total Vehicle Mass, S-IVB First Burn Phase, Kilograms.
  Mass total_vehicle_at_s_ivb_1st_90_percent_thrust = 161143 * Kilogram;
  Mass total_vehicle_at_s_ivb_1st_eco               = 128095 * Kilogram;

  // Table 21-7. Total Vehicle Mass, S-IVB Second Burn Phase, Kilograms.
  Mass total_vehicle_at_s_ivb_2nd_90_percent_thrust = 126780 * Kilogram;
  Mass total_vehicle_at_s_ivb_2nd_eco               =  59285 * Kilogram;

  Vector<double, World> e_y({0, 1, 0});
  Manœuvre<World> first_burn(thrust_1st,
                             total_vehicle_at_s_ivb_1st_90_percent_thrust,
                             specific_impulse_1st, e_y);
  EXPECT_EQ(lox_flowrate_1st + fuel_flowrate_1st, first_burn.mass_flow());

  first_burn.set_duration(s_ivb_1st_eco - s_ivb_1st_90_percent_thrust);
  EXPECT_EQ(total_vehicle_at_s_ivb_1st_eco, first_burn.final_mass());
  EXPECT_EQ(1 * Metre / Second, first_burn.Δv());

  first_burn.set_initial_time(s_ivb_1st_90_percent_thrust);
  EXPECT_EQ(s_ivb_1st_eco, first_burn.final_time());
  EXPECT_THAT(
      first_burn.acceleration()(first_burn.initial_time()).Norm(),
      AllOf(Gt(5 * Metre / Pow<2>(Second)), Lt(6.25 * Metre / Pow<2>(Second))));
  EXPECT_THAT(first_burn.acceleration()(first_burn.initial_time()).Norm(),
              AllOf(Gt(6.25 * Metre / Pow<2>(Second)),
                    Lt(7.5 * Metre / Pow<2>(Second))));
  EXPECT_THAT(first_burn.acceleration()(range_zero + 600 * Second).Norm(),
              AllOf(Gt(6.25 * Metre / Pow<2>(Second)),
                    Lt(6.25 * Metre / Pow<2>(Second))));
}

}  // namespace ksp_plugin
}  // namespace principia
