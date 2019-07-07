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
using ::testing::Property;
using ::testing::Lt;

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
  AngularFrequency const Ωʹ_tp = -2.076659 * Degree / Day;
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  112.43 * Minute, Ωʹ_tp, *earth_, /*max_abs_Cᴛₒ=*/50),
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
  OrbitRecurrence const метеор_3_7_recurrence(13, +7, 71);
  EXPECT_THAT(
      AbsoluteError(-27.48 * Degree, метеор_3_7_recurrence.equatorial_shift()),
      Lt(0.01 * Degree));
}

}  // namespace astronomy
}  // namespace principia
