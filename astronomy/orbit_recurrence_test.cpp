#include "astronomy/orbit_recurrence.hpp"

#include <limits>
#include <memory>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {

using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Lt;
using ::testing::Property;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_orbit_recurrence;
using namespace principia::base::_not_null;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics;

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

// This test exercises the Floating–integral conversions in cases where the
// behaviour is undefined.
TEST_F(OrbitRecurrenceTest, NearbyInt) {
  {
    int const i = std::nearbyint(0.1);
    EXPECT_EQ(0, i);
  }
#if PRINCIPIA_COMPILER_MSVC
  {
    int const i = std::nearbyint(1.1e20);
    EXPECT_EQ(-2147483648, i);
  }
  {
    int const i = std::nearbyint(std::numeric_limits<double>::quiet_NaN());
    EXPECT_EQ(-2147483648, i);
  }
  {
    int const i = std::nearbyint(std::numeric_limits<double>::infinity());
    EXPECT_EQ(-2147483648, i);
  }
#endif
  {
    using namespace principia::astronomy::_orbit_recurrence::internal;
    auto const status_or_i = SafeNearbyInt(0.1);
    EXPECT_THAT(status_or_i, IsOkAndHolds(0));
  }
  {
    using namespace principia::astronomy::_orbit_recurrence::internal;
    auto const status_or_i = SafeNearbyInt(1.1e20);
    EXPECT_THAT(status_or_i.status(),
                StatusIs(absl::StatusCode::kInvalidArgument));
  }
  {
    using namespace principia::astronomy::_orbit_recurrence::internal;
    auto const status_or_i =
        SafeNearbyInt(std::numeric_limits<double>::quiet_NaN());
    EXPECT_THAT(status_or_i.status(),
                StatusIs(absl::StatusCode::kInvalidArgument));
  }
  {
    using namespace principia::astronomy::_orbit_recurrence::internal;
    auto const status_or_i =
        SafeNearbyInt(std::numeric_limits<double>::infinity());
    EXPECT_THAT(status_or_i.status(),
                StatusIs(absl::StatusCode::kInvalidArgument));
  }
}

TEST_F(OrbitRecurrenceTest, ClosestRecurrence) {
  // Orbits from example 11.13.
  // The values for the nodal precession and for the period of Z-Earth are from
  // http://climserv.ipsl.polytechnique.fr/ixion/ (the value given for the
  // period of Z-Earth in figure 11.16(c) has insufficient precision).
  AngularFrequency const Ωʹꜱ = 0.985647 * Degree / Day;
  // SPOT-4, figure 11.15(a).
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  101.46 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/50),
              IsOkAndHolds(AllOf(Property(&OrbitRecurrence::νₒ, 14),
                                 Property(&OrbitRecurrence::Dᴛₒ, 5),
                                 Property(&OrbitRecurrence::Cᴛₒ, 26))));
  // Terra, figure 11.15(b).
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  98.88 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/50),
              IsOkAndHolds(AllOf(Property(&OrbitRecurrence::νₒ, 15),
                                 Property(&OrbitRecurrence::Dᴛₒ, -7),
                                 Property(&OrbitRecurrence::Cᴛₒ, 16))));

  // TOPEX/Poséidon, figure 11.16(b).
  AngularFrequency const Ωʹ_topex_poséidon = -2.076659 * Degree / Day;
  EXPECT_THAT(
      OrbitRecurrence::ClosestRecurrence(
          112.43 * Minute, Ωʹ_topex_poséidon, *earth_, /*max_abs_Cᴛₒ=*/50),
      IsOkAndHolds(AllOf(Property(&OrbitRecurrence::νₒ, 13),
                         Property(&OrbitRecurrence::Dᴛₒ, -3),
                         Property(&OrbitRecurrence::Cᴛₒ, 10))));

  // Z-Earth, figure 11.16(c).
  EXPECT_THAT(OrbitRecurrence::ClosestRecurrence(
                  95.097610 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/300),
              IsOkAndHolds(AllOf(Property(&OrbitRecurrence::νₒ, 15),
                                 Property(&OrbitRecurrence::Dᴛₒ, 39),
                                 Property(&OrbitRecurrence::Cᴛₒ, 274))));
  // Limiting Cᴛₒ to 50 days, we find the 7-day subcycle instead.
  EXPECT_THAT(
      OrbitRecurrence::ClosestRecurrence(
          95.097610 * Minute, Ωʹꜱ, *earth_, /*max_abs_Cᴛₒ=*/50),
      IsOkAndHolds(AllOf(Property(&OrbitRecurrence::number_of_revolutions, 106),
                         Property(&OrbitRecurrence::Cᴛₒ, 7))));
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

TEST_F(OrbitRecurrenceTest, RetrogradeRotation) {
  // Capderou does not study ground track recurrence around these bodies; the
  // expected values for the recurrence triple are nothing more than a change
  // detector.
  // Example 16.1.
  auto const status_or_magellan =
      OrbitRecurrence::ClosestRecurrence(91.914680 * Minute,
                                         -0.000464 * Degree / Day,
                                         *venus_,
                                         /*max_abs_Cᴛₒ=*/50);
  ASSERT_OK(status_or_magellan);
  auto const& magellan = status_or_magellan.value();
  EXPECT_THAT(AbsoluteError(0.094526 * Degree, magellan.equatorial_shift()),
              Lt(0.000001 * Degree));
  // There are approximately 3807 orbits per sidereal day; this value of νₒ is
  // consistent with the low precession.
  EXPECT_THAT(magellan,
              AllOf(Property(&OrbitRecurrence::νₒ, -3808),
                    Property(&OrbitRecurrence::Dᴛₒ, 1),
                    Property(&OrbitRecurrence::Cᴛₒ, -2)));
  EXPECT_THAT(magellan.subcycle(), Eq(-1));
  // Example 16.9.
  auto const status_or_triton_orbiter =
      OrbitRecurrence::ClosestRecurrence(169.584671 * Minute,
                                         -0.074331 * Degree / Day,
                                         *triton_,
                                         /*max_abs_Cᴛₒ=*/50);
  ASSERT_OK(status_or_triton_orbiter);
  auto const& triton_orbiter = status_or_triton_orbiter.value();
  EXPECT_THAT(AbsoluteError(7.2 * Degree, triton_orbiter.equatorial_shift()),
              Lt(0.1 * Degree));
  EXPECT_THAT(triton_orbiter,
              AllOf(Property(&OrbitRecurrence::νₒ, -50),
                    Property(&OrbitRecurrence::Dᴛₒ, -1),
                    Property(&OrbitRecurrence::Cᴛₒ, -27)));
  EXPECT_THAT(triton_orbiter.subcycle(), Eq(-1));
}

TEST_F(OrbitRecurrenceTest, OneDayRecurrenceCycle) {
  AngularFrequency const Ωʹꜱ = 0.985647 * Degree / Day;
  // Example 11.10: FORMOSAT-2, formerly ROCSAT-2.
  auto const status_or_福爾摩沙衛星二號 =
      OrbitRecurrence::ClosestRecurrence(102.86 * Minute,
                                         Ωʹꜱ,
                                         *earth_,
                                         /*max_abs_Cᴛₒ=*/50);
  ASSERT_OK(status_or_福爾摩沙衛星二號);
  auto const& 福爾摩沙衛星二號 = status_or_福爾摩沙衛星二號.value();
  EXPECT_THAT(福爾摩沙衛星二號,
              AllOf(Property(&OrbitRecurrence::νₒ, 14),
                    Property(&OrbitRecurrence::Dᴛₒ, 0),
                    Property(&OrbitRecurrence::Cᴛₒ, 1)));
  EXPECT_THAT(福爾摩沙衛星二號.grid_interval(),
              Eq(福爾摩沙衛星二號.base_interval()));
  EXPECT_THAT(福爾摩沙衛星二號.subcycle(), Eq(0));
}

}  // namespace astronomy
}  // namespace principia
