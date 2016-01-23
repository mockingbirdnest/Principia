
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
using geometry::JulianDate;
using geometry::Sign;
using integrators::McLachlanAtela1992Order5Optimal;
using numerics::Bisect;
using physics::Ephemeris;
using quantities::ArcCos;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using testing_utilities::AbsoluteError;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

namespace {

Time const bisection_interval = 10 * Minute;

Sign const U14 = Sign(1);
Sign const U23 = Sign(-1);

// MassiveBody will need radius information.  Or non-hardcoded data from
// https://github.com/mockingbirdnest/Principia/blob/master/astronomy/gravity_model.proto.txt

// Interesting analysis of the Sun radius at http://www.icra.it/solar/Sole2.pdf.
Length const r_sun = 696000.0 * Kilo(Metre);
Length const r_earth = 6378.1363 * Kilo(Metre);
Length const r_moon = 1738.0 * Kilo(Metre);

}  // namespace

class EclipseTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
        45 * Minute, 5 * Milli(Metre));
    atmospheric_depth_ = ComputeAtmosphericDepthAtFirstContact();
  }

  static Length ComputeAtmosphericDepthAtFirstContact() {
    SolarSystem<ICRFJ2000Equator> solar_system_first_contact;
    solar_system_first_contact.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2451564_587154910.proto.txt");

    auto const q_sun_first_contact =
        solar_system_first_contact.initial_state("Sun").position();
    auto const q_earth_first_contact =
        solar_system_first_contact.initial_state("Earth").position();
    auto const q_moon_first_contact =
        solar_system_first_contact.initial_state("Moon").position();

    // Angle between the Sun-Earth axis and the tangent ray from the Sun to the
    // Earth.
    auto const alpha = [&q_earth_first_contact,
                        &q_sun_first_contact](Length const& depth) {
      return ArcSin((r_sun + r_earth + depth) /
                    (q_sun_first_contact - q_earth_first_contact).Norm());
    };

    // Angle between the Earth-Moon axis and the tangent ray from the Earth to
    // the moon.
    auto const beta = [&q_earth_first_contact,
                       &q_moon_first_contact](Length const& depth) {
      return ArcSin((r_moon + r_earth + depth) /
                    (q_moon_first_contact - q_earth_first_contact).Norm());
    };

    // Angle between the Sun-Earth axis and the Earth-Moon axis.
    auto const gamma =
        ArcCos(InnerProduct(q_sun_first_contact - q_earth_first_contact,
                            q_earth_first_contact - q_moon_first_contact) /
               ((q_earth_first_contact - q_moon_first_contact).Norm() *
                (q_sun_first_contact - q_earth_first_contact).Norm()));

    // Find the atmospheric depth that cancels the error at this contact.
    Length const& actual_depth = Bisect(
        [alpha, beta, gamma](Length const& depth) {
          return gamma - alpha(depth) - beta(depth);
        },
        0 * Kilo(Metre),
        1000 * Kilo(Metre));
    return actual_depth;
  }

  // A positive |time_error| means that the actual contact happens after
  // |current_time|.
  void CheckLunarUmbralEclipse(Instant const& current_time,
                               Sign const moon_offset_sign,
                               Angle const& angular_error,
                               Time const& time_error) {
    ephemeris_->Prolong(current_time + bisection_interval);
    auto const sun = solar_system_1950_.massive_body(*ephemeris_, "Sun");
    auto const earth = solar_system_1950_.massive_body(*ephemeris_, "Earth");
    auto const moon = solar_system_1950_.massive_body(*ephemeris_, "Moon");

    // Angle formed by a right circular cone with sides defined by tangent lines
    // between Sun and Earth, and axis running through the centers of each.
    auto const umbral_half_aperture = [this, earth, moon, sun](
        Instant const& t) {
      auto const q_sun =
          ephemeris_->trajectory(sun)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_moon =
          ephemeris_->trajectory(moon)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_earth =
          ephemeris_->trajectory(earth)->EvaluatePosition(t, /*hint=*/nullptr);
      return ArcSin((r_sun - (r_earth + atmospheric_depth_)) /
                    (q_sun - q_earth).Norm());
    };

    auto const earth_moon_angle = [this,
                                   earth,
                                   moon,
                                   moon_offset_sign,
                                   sun,
                                   umbral_half_aperture](Instant const& t) {
      auto const q_sun =
          ephemeris_->trajectory(sun)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_moon =
          ephemeris_->trajectory(moon)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_earth =
          ephemeris_->trajectory(earth)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const apex_of_moon_locus_at_umbral_contact =
          q_earth +
          Normalize(q_earth - q_sun) *
              (r_earth + atmospheric_depth_ + moon_offset_sign * r_moon) /
              Sin(umbral_half_aperture(t));
      // Angle between Earth and Moon as seen at
      // |apex_of_moon_locus_at_umbral_contact|.
      return ArcCos(
          InnerProduct(apex_of_moon_locus_at_umbral_contact - q_earth,
                       apex_of_moon_locus_at_umbral_contact - q_moon) /
          ((apex_of_moon_locus_at_umbral_contact - q_moon).Norm() *
           (apex_of_moon_locus_at_umbral_contact - q_earth).Norm()));
    };

    // We are at the desired contact if the angle between Earth and Moon from
    // the apex of locus of the moon at that contact is the same value as the
    // half-aperture of the umbra (Earth-Sun cone).
    EXPECT_THAT(AbsoluteError(umbral_half_aperture(current_time),
                              earth_moon_angle(current_time)),
                AllOf(Lt(angular_error), Gt(0.5 * angular_error)))
        << NAMED(umbral_half_aperture(current_time)) << ", "
        << NAMED(earth_moon_angle(current_time)) << ", " << NAMED(current_time);

    Instant const& actual_contact_time = Bisect(
        [earth_moon_angle, umbral_half_aperture](Instant const& t) {
          return umbral_half_aperture(t) - earth_moon_angle(t);
        },
        current_time - bisection_interval,
        current_time + bisection_interval);
    EXPECT_EQ(Sign(actual_contact_time - current_time),
              Sign(time_error)) << NAMED(time_error);
    EXPECT_THAT(AbsoluteError(actual_contact_time, current_time),
                AllOf(Lt(Abs(time_error)), Gt(0.5 * Abs(time_error))))
        << NAMED(actual_contact_time) << ", " << NAMED(current_time);
  }

  // A positive |time_error| means that the actual contact happens after
  // |current_time|.
  void CheckLunarPenumbralEclipse(Instant const& current_time,
                                  Sign const moon_offset_sign,
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
      return ArcSin((r_sun + r_earth + atmospheric_depth_) /
                    (q_sun - q_earth).Norm());
    };

    auto const earth_moon_angle = [this,
                                   earth,
                                   moon,
                                   moon_offset_sign,
                                   sun,
                                   penumbral_half_aperture](Instant const& t) {
      auto const q_sun =
          ephemeris_->trajectory(sun)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_moon =
          ephemeris_->trajectory(moon)->EvaluatePosition(t, /*hint=*/nullptr);
      auto const q_earth =
          ephemeris_->trajectory(earth)->EvaluatePosition(t, /*hint=*/nullptr);

      auto const apex_of_moon_locus_at_penumbral_contact =
          q_earth +
          Normalize(q_sun - q_earth) *
              (r_earth + atmospheric_depth_ + moon_offset_sign * r_moon) /
              Sin(penumbral_half_aperture(t));
      // Angle between Earth and Moon as seen at
      // apex_of_moon_locus_at_penumbral_contact.
      return ArcCos(
          InnerProduct(apex_of_moon_locus_at_penumbral_contact - q_earth,
                       apex_of_moon_locus_at_penumbral_contact - q_moon) /
          ((apex_of_moon_locus_at_penumbral_contact - q_moon).Norm() *
           (apex_of_moon_locus_at_penumbral_contact - q_earth).Norm()));
    };

    // We are at the desired contact if the angle between Earth and Moon from
    // the apex of locus of the moon at that contact is the same value as the
    // half-aperture of the penumbra.
    EXPECT_THAT(AbsoluteError(penumbral_half_aperture(current_time),
                              earth_moon_angle(current_time)),
                AllOf(Lt(angular_error), Gt(0.5 * angular_error)))
        << NAMED(penumbral_half_aperture(current_time)) << ", "
        << NAMED(earth_moon_angle(current_time)) << ", " << NAMED(current_time);

    Instant const& actual_contact_time = Bisect(
        [earth_moon_angle, penumbral_half_aperture](Instant const& t) {
          return penumbral_half_aperture(t) - earth_moon_angle(t);
        },
        current_time - bisection_interval,
        current_time + bisection_interval);
    EXPECT_EQ(Sign(actual_contact_time - current_time),
              Sign(time_error)) << NAMED(time_error);
    EXPECT_THAT(AbsoluteError(actual_contact_time, current_time),
                AllOf(Lt(Abs(time_error)), Gt(0.5 * Abs(time_error))))
        << NAMED(actual_contact_time) << ", " << NAMED(current_time);
  }

  static SolarSystem<ICRFJ2000Equator> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  static Length atmospheric_depth_;
};

SolarSystem<ICRFJ2000Equator> EclipseTest::solar_system_1950_;
std::unique_ptr<Ephemeris<ICRFJ2000Equator>> EclipseTest::ephemeris_;
Length EclipseTest::atmospheric_depth_;

#if !defined(_DEBUG)

TEST_F(EclipseTest, AtmosphericDepth) {
  EXPECT_THAT(atmospheric_depth_, AllOf(Lt(92 * Kilo(Metre)),
                                        Gt(91 * Kilo(Metre))));
}

TEST_F(EclipseTest, Year1950) {
  // Times are TDB Julian Day for 1950-04-02.
  auto P1 = JulianDate(2433374.25788409);  // 18:10:49 UT
  auto U1 = JulianDate(2433374.29850909);  // 19:09:19
  auto U2 = JulianDate(2433374.354979);    // 20:30:38
  auto U3 = JulianDate(2433374.37367113);  // 20:57:33
  auto U4 = JulianDate(2433374.43016419);  // 22:18:54
  auto P4 = JulianDate(2433374.47075446);  // 23:17:21

  CheckLunarPenumbralEclipse(P1, U14, 5E-7 * Radian,   -1 * Second);
  CheckLunarUmbralEclipse(U1, U14,    5E-7 * Radian,   -1 * Second);
  CheckLunarUmbralEclipse(U2, U23,    3E-5 * Radian,  -60 * Second);
  CheckLunarUmbralEclipse(U3, U23,    5E-5 * Radian,  120 * Second);
  CheckLunarUmbralEclipse(U4, U14,    4E-5 * Radian,   60 * Second);
  CheckLunarPenumbralEclipse(P4, U14, 3E-5 * Radian,   60 * Second);

  // Times are TDB Julian Day for 1950-09-26.
  P1 = JulianDate(2433550.55712016);  // 01:21:43 UT
  U1 = JulianDate(2433550.60578913);  // 02:31:48
  U2 = JulianDate(2433550.66325441);  // 03:54:33
  U3 = JulianDate(2433550.69399515);  // 04:38:49
  U4 = JulianDate(2433550.75144885);  // 06:01:33
  P4 = JulianDate(2433550.800222);    // 07:11:47

  CheckLunarPenumbralEclipse(P1, U14, 2E-6 * Radian,    3 * Second);
  CheckLunarUmbralEclipse(U1, U14,    2E-6 * Radian,    4 * Second);
  CheckLunarUmbralEclipse(U2, U23,    2E-5 * Radian,  -30 * Second);
  CheckLunarUmbralEclipse(U3, U23,    6E-5 * Radian,  100 * Second);
  CheckLunarUmbralEclipse(U4, U14,    5E-5 * Radian,   70 * Second);
  CheckLunarPenumbralEclipse(P4, U14, 3E-5 * Radian,   70 * Second);
}

TEST_F(EclipseTest, Year1951) {
  // Times are TDB Julian Day for 1951-03-23.
  auto P1 = JulianDate(2433728.86842806);  // 08:50:50
  auto P4 = JulianDate(2433729.01725909);  // 12:24:19

  CheckLunarPenumbralEclipse(P1, U14, 5E-6 * Radian, -15 * Second);
  CheckLunarPenumbralEclipse(P4, U14, 3E-5 * Radian,  70 * Second);

  // Times are TDB Julian Day for 1951-09-15.
  P1 = JulianDate(2433904.93736321);  // 10:29:16
  P4 = JulianDate(2433905.1002799);   // 14:23:52

  CheckLunarPenumbralEclipse(P1, U14, 5E-6 * Radian, -13 * Second);
  CheckLunarPenumbralEclipse(P4, U14, 3E-5 * Radian,  70 * Second);
}

TEST_F(EclipseTest, Year1952) {
  // Times are TDB Julian Day for 1952-02-11 (or 10 for P1).
  auto P1 = JulianDate(2434053.42282623);  // P1 = 22:08:20 UT
  auto U1 = JulianDate(2434053.50334705);  // U1 = 00:04:17
  auto U4 = JulianDate(2434053.55203917);  // U4 = 01:14:24
  auto P4 = JulianDate(2434053.63249055);  // P4 = 03:10:15

  CheckLunarPenumbralEclipse(P1, U14, 3E-6 * Radian,  -9 * Second);
  CheckLunarUmbralEclipse(U1, U14,    2E-5 * Radian, -90 * Second);
  CheckLunarUmbralEclipse(U4, U14,    4E-5 * Radian, 150 * Second);
  CheckLunarPenumbralEclipse(P4, U14, 3E-5 * Radian,  80 * Second);

  // Times are TDB Julian Day for 1952-08-05.
  P1 = JulianDate(2434230.22830075);  // P1 = 17:28:13 UT
  U1 = JulianDate(2434230.27385631);  // U1 = 18:33:49
  U4 = JulianDate(2434230.37606695);  // U4 = 21:01:00
  P4 = JulianDate(2434230.42161093);  // P4 = 22:06:35

  CheckLunarPenumbralEclipse(P1, U14, 6E-6 * Radian, -13 * Second);
  CheckLunarUmbralEclipse(U1, U14,    2E-5 * Radian, -30 * Second);
  CheckLunarUmbralEclipse(U4, U14,    4E-5 * Radian,  70 * Second);
  CheckLunarPenumbralEclipse(P4, U14, 3E-5 * Radian,  60 * Second);

  // Later on for additional accuracy: 2 * ArcTan((x_norm_y -
  // y_normx).Norm(),(x_norm_y + y_norm_x).Norm())
  // x_norm_y = x * y.Norm() and y_norm_x = y * x.Norm()
}

#if 0
TEST_F(EclipseTest, Year2048) {
  // Times are TDB Julian Day for 2048-01-01.
  auto P1 = JulianDate(2469076.66235167);  // 03:52:39 UT
  auto U1 = JulianDate(2469076.71279148);  // 05:05:17
  auto U2 = JulianDate(2469076.76776833);  // 06:24:27
  auto U3 = JulianDate(2469076.80661092);  // 07:20:23
  auto U4 = JulianDate(2469076.86158778);  // 08:39:33
  auto P4 = JulianDate(2469076.91195815);  // 09:52:05

  CheckLunarPenumbralEclipse(P1, U14, 2E-4 * Radian, -370 * Second);
  CheckLunarUmbralEclipse(U1, U14,    3E-4 * Radian, -370 * Second);
  CheckLunarUmbralEclipse(U2, U23,    3E-4 * Radian, -390 * Second);
  CheckLunarUmbralEclipse(U3, U23,    3E-4 * Radian, -290 * Second);
  CheckLunarUmbralEclipse(U4, U14,    2E-4 * Radian, -310 * Second);
  CheckLunarPenumbralEclipse(P4, U14, 2E-4 * Radian, -310 * Second);
}
#endif

#endif

}  // namespace physics
}  // namespace principia
