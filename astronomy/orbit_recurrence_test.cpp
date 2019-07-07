#include "astronomy/orbit_recurrence.hpp"

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/solar_system.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {

using base::not_null;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::AngularFrequency;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Minute;
using testing_utilities::AbsoluteError;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Lt;
using ::testing::Property;

class OrbitRecurrenceTest : public ::testing::Test {
 protected:
  OrbitRecurrenceTest()
      : solar_system_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt"),
        earth_(SolarSystem<ICRS>::MakeRotatingBody(
            solar_system_.gravity_model_message("Earth"))),
        venus_(SolarSystem<ICRS>::MakeRotatingBody(
            solar_system_.gravity_model_message("Venus"))),
        triton_(SolarSystem<ICRS>::MakeRotatingBody(
            solar_system_.gravity_model_message("Triton"))) {}

  SolarSystem<ICRS> solar_system_;
  not_null<std::unique_ptr<RotatingBody<ICRS>>> earth_;
  not_null<std::unique_ptr<RotatingBody<ICRS>>> venus_;
  not_null<std::unique_ptr<RotatingBody<ICRS>>> triton_;
};

TEST_F(OrbitRecurrenceTest, ClosestRecurrence) {
  // Orbits from example 11.13.
  // The values for the nodal precession and for the period of Z-Earth are from
  // http://climserv.ipsl.polytechnique.fr/ixion/ (the value given for the
  // period of Z-Earth in figure 11.16(c) has insufficient precision).
  AngularFrequency const Ωʹꜱ = 0.985647 * Degree / Day;
  // SPOT-4, figure 11.15(a).
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  101.46 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/50),
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 5),
                    Property(&OrbitRecurrence::Cᴛₒ, 26)));
  // Terra, figure 11.15(b).
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  98.88 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/50),
              AllOf(Property(&OrbitRecurrence::νₒ, 15),
                    Property(&OrbitRecurrence::Dᴛₒ, -7),
                    Property(&OrbitRecurrence::Cᴛₒ, 16)));

  // TOPEX/Poséidon, figure 11.16(b).
  AngularFrequency const Ωʹ_topex_poséidon = -2.076659 * Degree / Day;
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
          112.43 * Minute, Ωʹ_topex_poséidon, *earth_, /*max_abs_Cᴛₒ=*/50),
      AllOf(Property(&OrbitRecurrence::νₒ, 13),
            Property(&OrbitRecurrence::Dᴛₒ, -3),
            Property(&OrbitRecurrence::Cᴛₒ, 10)));

  // Z-Earth, figure 11.16(c).
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  95.097610 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/300),
              AllOf(Property(&OrbitRecurrence::νₒ, 15),
                    Property(&OrbitRecurrence::Dᴛₒ, 39),
                    Property(&OrbitRecurrence::Cᴛₒ, 274)));
  // Limiting Cᴛₒ to 50 days, we find the 7-day subcycle instead.
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  95.097610 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/50),
              AllOf(Property(&OrbitRecurrence::number_of_revolutions, 106),
                    Property(&OrbitRecurrence::Cᴛₒ, 7)));
}

TEST_F(OrbitRecurrenceTest, EquatorialShift) {
  // Example 8.2: Метеор-3 № 7.
  OrbitRecurrence const метеор_3_7(13, +7, 71);
  EXPECT_THAT(AbsoluteError(-27.48 * Degree, метеор_3_7.equatorial_shift()),
      Lt(0.01 * Degree));
}

TEST_F(OrbitRecurrenceTest, GridInterval) {
  // Example 11.8: TOPEX/Poséidon.
  OrbitRecurrence const topex_poséidon(13, -3, 10);
  EXPECT_THAT(AbsoluteError(2.8346 * Degree, topex_poséidon.grid_interval()),
              Lt(0.0001 * Degree));
  EXPECT_THAT(AbsoluteError(28.35 * Degree, topex_poséidon.base_interval()),
              Lt(0.01 * Degree));
}

TEST_F(OrbitRecurrenceTest, Subcycle) {
  // Example 11.9.
  OrbitRecurrence const spot_5(14, +5, 26);
  EXPECT_THAT(spot_5.subcycle(), Eq(5));
  OrbitRecurrence const terra(15, -7, 16);
  EXPECT_THAT(terra.subcycle(), Eq(7));
  OrbitRecurrence const adeos_1(14, +11, 41);
  EXPECT_THAT(adeos_1.subcycle(), Eq(15));
  OrbitRecurrence const topex_poséidon(13, -3, 10);
  EXPECT_THAT(topex_poséidon.subcycle(), Eq(3));
  OrbitRecurrence const icesat(15, -22, 183);
  EXPECT_THAT(icesat.subcycle(), Eq(25));
}

}  // namespace astronomy
}  // namespace principia
