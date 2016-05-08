
#include "geometry/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "numerics/root_finders.hpp"
#include "physics/solar_system.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::AngleBetween;
using geometry::JulianDate;
using geometry::Sign;
using integrators::McLachlanAtela1992Order5Optimal;
using numerics::Bisect;
using physics::Ephemeris;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using testing_utilities::AbsoluteError;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

namespace {

Time const bisection_interval = 10 * Minute;

Sign const U14 = Sign(1);
Sign const U23 = Sign(-1);

char const arrow[] = "-------------------> ";

}  // namespace

class EclipseTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
            McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
            /*step=*/45 * Minute));
    r_sun_ = solar_system_1950_.mean_radius("Sun");
    r_moon_ = solar_system_1950_.mean_radius("Moon");
    r_earth_ = solar_system_1950_.mean_radius("Earth");
  }

  // A positive |time_error| means that the actual contact happens after
  // |current_time|.
  void CheckSolarUmbralEclipse(Instant const& current_time,
                               Sign const earth_offset_sign,
                               Angle const& angular_error,
                               Time const& time_error) {
    ephemeris_->Prolong(current_time + bisection_interval);
    auto const sun = solar_system_1950_.massive_body(*ephemeris_, "Sun");
    auto const earth = solar_system_1950_.massive_body(*ephemeris_, "Earth");
    auto const moon = solar_system_1950_.massive_body(*ephemeris_, "Moon");

    // Angle formed by a right circular cone with sides defined by tangent lines
    // between Sun and Moon, and axis running through the centers of each.
    auto const umbral_half_aperture = [this, earth, moon, sun](
        Instant const& t) {
      auto const q_sun =
          ephemeris_->trajectory(sun)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_moon =
          ephemeris_->trajectory(moon)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_earth =
          ephemeris_->trajectory(earth)->EvaluatePosition(t, /*hint=*/nullptr);
      return ArcSin((r_sun_ - r_moon_) /
                    (q_sun - q_moon).Norm());
    };

    auto const moon_earth_angle = [this,
                                   earth,
                                   moon,
                                   earth_offset_sign,
                                   sun,
                                   umbral_half_aperture](Instant const& t) {
      auto const q_sun =
          ephemeris_->trajectory(sun)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_moon =
          ephemeris_->trajectory(moon)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_earth =
          ephemeris_->trajectory(earth)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const apex_of_earth_locus_at_umbral_contact =
          q_moon +
          Normalize(q_moon - q_sun) *
              (r_moon_ + earth_offset_sign * r_earth_) /
              Sin(umbral_half_aperture(t));
      // Angle between Earth and Moon as seen at
      // |apex_of earth_locus_at_umbral_contact|.
      return AngleBetween(apex_of earth_locus_at_umbral_contact - q_moon,
                          apex_of earth_locus_at_umbral_contact - q_earth);
    };

    // We are at the desired contact if the angle between Earth and Moon from
    // the apex of locus of the moon at that contact is the same value as the
    // half-aperture of the umbra (Earth-Sun cone).
    EXPECT_THAT(AbsoluteError(umbral_half_aperture(current_time),
                              moon_earth_angle(current_time)),
                AllOf(Lt(angular_error), Gt(0.5 * angular_error)))
        << NAMED(umbral_half_aperture(current_time)) << ", "
        << NAMED(moon_earth_angle(current_time)) << ", " << NAMED(current_time);

    Instant const& actual_contact_time = Bisect(
        [moon_earth_angle, umbral_half_aperture](Instant const& t) {
          return umbral_half_aperture(t) - moon_earth_angle(t);
        },
        current_time - bisection_interval,
        current_time + bisection_interval);
    EXPECT_THAT(Sign(actual_contact_time - current_time), Eq(Sign(time_error)))
        << NAMED(actual_contact_time - current_time);
    EXPECT_THAT(AbsoluteError(actual_contact_time, current_time),
                AllOf(Lt(Abs(time_error)), Gt(0.9 * Abs(time_error))))
        << NAMED(actual_contact_time) << ", " << NAMED(current_time);
    LOG(INFO) << arrow << AbsoluteError(umbral_half_aperture(current_time),
                                        moon_earth_angle(current_time))
              << " " << actual_contact_time - current_time;
  }

  // A positive |time_error| means that the actual contact happens after
  // |current_time|.
  void CheckSolarPenumbralEclipse(Instant const& current_time,
                                  Sign const earth_offset_sign,
                                  Angle const& angular_error,
                                  Time const& time_error) {
    ephemeris_->Prolong(current_time + bisection_interval);
    auto const sun = solar_system_1950_.massive_body(*ephemeris_, "Sun");
    auto const earth = solar_system_1950_.massive_body(*ephemeris_, "Earth");
    auto const moon = solar_system_1950_.massive_body(*ephemeris_, "Moon");

    auto const penumbral_half_aperture = [this, earth, moon, sun](
        Instant const& t) {
      auto const q_sun =
          ephemeris_->trajectory(sun)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_moon =
          ephemeris_->trajectory(moon)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_earth =
          ephemeris_->trajectory(earth)->EvaluatePosition(t, /*hint=*/nullptr);
      return ArcSin((r_sun_ + r_moon_) / (q_sun - q_moon).Norm());
    };

    auto const moon_earth_angle = [this,
                                   earth,
                                   moon,
                                   earth_offset_sign,
                                   sun,
                                   penumbral_half_aperture](Instant const& t) {
      auto const q_sun =
          ephemeris_->trajectory(sun)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_moon =
          ephemeris_->trajectory(moon)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_earth =
          ephemeris_->trajectory(earth)->EvaluatePosition(t, /*hint=*/nullptr);

      auto const apex_of_earth_locus_at_penumbral_contact =
          q_moon +
          Normalize(q_sun - moon) *
              (r_moon_ + earth_offset_sign * r_earth_) /
              Sin(penumbral_half_aperture(t));
      // Angle between Moon and Earth as seen at
      // apex_of_earth_locus_at_penumbral_contact.
      return AngleBetween(apex_of_earth_locus_at_penumbral_contact - q_moon,
                          apex_of_earth_locus_at_penumbral_contact - q_earth);
    };

    // We are at the desired contact if the angle between Moon and Earth from
    // the apex of locus of the earth at that contact is the same value as the
    // half-aperture of the penumbra.
    EXPECT_THAT(AbsoluteError(penumbral_half_aperture(current_time),
                              moon_earth_angle(current_time)),
                AllOf(Lt(angular_error), Gt(0.5 * angular_error)))
        << NAMED(penumbral_half_aperture(current_time)) << ", "
        << NAMED(moon_earth_angle(current_time)) << ", " << NAMED(current_time);

    Instant const& actual_contact_time = Bisect(
        [moon_earth_angle, penumbral_half_aperture](Instant const& t) {
          return penumbral_half_aperture(t) - moon_earth_angle(t);
        },
        current_time - bisection_interval,
        current_time + bisection_interval);
    EXPECT_THAT(Sign(actual_contact_time - current_time), Eq(Sign(time_error)))
        << NAMED(actual_contact_time - current_time);
    EXPECT_THAT(AbsoluteError(actual_contact_time, current_time),
                AllOf(Lt(Abs(time_error)), Gt(0.9 * Abs(time_error))))
        << NAMED(actual_contact_time) << ", " << NAMED(current_time);
    LOG(INFO) << arrow
              << Sign(penumbral_half_aperture(current_time) -
                      moon_earth_angle(current_time)) *
                     AbsoluteError(penumbral_half_aperture(current_time),
                                   moon_earth_angle(current_time))
              << " " << actual_contact_time - current_time;
  }

  static SolarSystem<ICRFJ2000Equator> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  static Length r_sun_;
  static Length r_earth_;
  static Length r_moon_;
};

SolarSystem<ICRFJ2000Equator> EclipseTest::solar_system_1950_;
std::unique_ptr<Ephemeris<ICRFJ2000Equator>> EclipseTest::ephemeris_;
Length EclipseTest::r_sun_;
Length EclipseTest::r_earth_;
Length EclipseTest::r_moon_;

#if !defined(_DEBUG)

TEST_F(EclipseTest, Year1950) {
  // Times are TDB Julian Day for 1950-03-18.
  auto P1 = JulianDate(2433374.04951604);  // 13:10:46.0 UT
  auto U1 = JulianDate(2433374.13129729);  // 15:08:31.9
  auto U4 = JulianDate(2433374.16370932);  // 15:55:12.3
  auto P4 = JulianDate(2433374.24546743);  // 17:52:56.2

  CheckSolarPenumbralEclipse(P1, U14, 2E-5 * Radian, 28 * Second);
  CheckSolarUmbralEclipse(U1, U14,    2E-5 * Radian, 30 * Second);
  CheckSolarUmbralEclipse(U4, U14,    2E-5 * Radian, 28 * Second);
  CheckSolarPenumbralEclipse(P4, U14, 2E-5 * Radian, 30 * Second);

  // Times are TDB Julian Day for 1950-09-12.
  P1 = JulianDate(2433536.55815026);  // 01:23:12.0 UT
  U1 = JulianDate(2433536.61806924);  // 02:49:29.0
  auto U2 = JulianDate(2433536.61896855);  // 02:50:46.7
  auto U3 = JulianDate(2433536.6852741);  // 04:26:15.5
  U4 = JulianDate(2433536.68622896  );  // 04:27:38.0
  P4 = JulianDate(2433536.7473204);  // 05:55:36.3

  CheckSolarPenumbralEclipse(P1, U14, 2E-5 * Radian, 34 * Second);
  CheckSolarUmbralEclipse(U1, U14,    3E-5 * Radian, 36 * Second);
  CheckSolarUmbralEclipse(U2, U23,    3E-5 * Radian, 42 * Second);
  CheckSolarUmbralEclipse(U3, U23,    2E-5 * Radian, 29 * Second);
  CheckSolarUmbralEclipse(U4, U14,    3E-5 * Radian, 34 * Second);
  CheckSolarPenumbralEclipse(P4, U14, 2E-5 * Radian, 36 * Second);
}

TEST_F(EclipseTest, Year1951) {
  // Times are TDB Julian Day for 1951-03-07.
  auto P1 = JulianDate(2433713.25311442);  // 18:03:56.9 UT
  auto U1 = JulianDate(2433713.29580192);  // 19:05:25.1
  auto U2 = JulianDate(2433713.29693617);  // 19:07:03.1
  auto U3 = JulianDate(2433713.44439335);  // 22:39:23.4
  auto U4 = JulianDate(2433713.44559242);  // 22:41:07.0 
  auto P4 = JulianDate(2433713.4882753);  // 23:42:34.8

  CheckSolarPenumbralEclipse(P1, U14, 9E-6 * Radian, 30 * Second);
  CheckSolarUmbralEclipse(U1, U14,    3E-5 * Radian, 36 * Second);
  CheckSolarUmbralEclipse(U2, U23,    3E-5 * Radian, 42 * Second);
  CheckSolarUmbralEclipse(U3, U23,    2E-5 * Radian, 29 * Second);
  CheckSolarUmbralEclipse(U4, U14,    3E-5 * Radian, 34 * Second);
  CheckSolarPenumbralEclipse(P4, U14, 8E-6 * Radian, 25 * Second);

  // Times are TDB Julian Day for 1951-09-01.
  P1 = JulianDate(2433890.91318846);  // P1 = 09:54:27.3 UT
  U1 = JulianDate(2433890.95685628);  // U1 = 10:57:20.2
  U2 = JulianDate(2433890.9587498);  // U2 = 11:00:03.8
  U3 = JulianDate(2433891.1133829);  // U3 = 14:42:44.1
  U4 = JulianDate(2433891.1152116);  // U4 = 14:45:22.1
  P4 = JulianDate(2433891.15882734);  // P4 = 15:48:10.5

  CheckSolarPenumbralEclipse(P1, U14, 9E-6 * Radian, 28 * Second);
  CheckSolarUmbralEclipse(U1, U14,    3E-5 * Radian, 36 * Second);
  CheckSolarUmbralEclipse(U2, U23,    3E-5 * Radian, 42 * Second);
  CheckSolarUmbralEclipse(U3, U23,    2E-5 * Radian, 29 * Second);
  CheckSolarUmbralEclipse(U4, U14,    3E-5 * Radian, 34 * Second);
  CheckSolarPenumbralEclipse(P4, U14, 8E-6 * Radian, 23 * Second);
}

#if 0
TEST_F(EclipseTest, Year2000) {
  // Times are TDB Julian Day for 2000-01-21.
  auto P1 = JulianDate(2451564.58715491);  // 02:04:26 UT
  auto U1 = JulianDate(2451564.62701602);  // 03:01:50
  auto U2 = JulianDate(2451564.67089334);  // 04:05:01
  auto U3 = JulianDate(2451564.72435399);  // 05:22:00
  auto U4 = JulianDate(2451564.76820815);  // 06:25:09
  auto P4 = JulianDate(2451564.80812714);  // 07:22:38

  CheckSolarPenumbralEclipse(P1, U14, 8E-5 * Radian, -167 * Second);
  CheckSolarUmbralEclipse(U1, U14,    2E-4 * Radian, -164 * Second);
  CheckSolarUmbralEclipse(U2, U23,    3E-4 * Radian, -167 * Second);
  CheckSolarUmbralEclipse(U3, U23,    2E-4 * Radian, -160 * Second);
  CheckSolarUmbralEclipse(U4, U14,    2E-4 * Radian, -161 * Second);
  CheckSolarPenumbralEclipse(P4, U14, 8E-5 * Radian, -160 * Second);
}

TEST_F(EclipseTest, Year2048) {
  // Times are TDB Julian Day for 2048-01-01.
  auto P1 = JulianDate(2469076.66235167);  // 03:52:39 UT
  auto U1 = JulianDate(2469076.71279148);  // 05:05:17
  auto U2 = JulianDate(2469076.76776833);  // 06:24:27
  auto U3 = JulianDate(2469076.80661092);  // 07:20:23
  auto U4 = JulianDate(2469076.86158778);  // 08:39:33
  auto P4 = JulianDate(2469076.91195815);  // 09:52:05

  CheckSolarPenumbralEclipse(P1, U14, 2E-4 * Radian, -338 * Second);
  CheckSolarUmbralEclipse(U1, U14,    3E-4 * Radian, -336 * Second);
  CheckSolarUmbralEclipse(U2, U23,    3E-4 * Radian, -335 * Second);
  CheckSolarUmbralEclipse(U3, U23,    3E-4 * Radian, -337 * Second);
  CheckSolarUmbralEclipse(U4, U14,    3E-4 * Radian, -336 * Second);
  CheckSolarPenumbralEclipse(P4, U14, 2E-4 * Radian, -335 * Second);
}
#endif

#endif

}  // namespace physics
}  // namespace principia
