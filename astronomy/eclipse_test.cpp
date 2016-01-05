
#include "geometry/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
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
using physics::Ephemeris;
using quantities::ArcCos;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Nano;
using testing_utilities::AbsoluteError;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

namespace {
Sign const U14 = Sign(1);
Sign const U23 = Sign(-1);
}  // namespace

class EclipseTest : public testing::Test {
 protected:
 /* EclipseTest() {
    solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
        45 * Minute, 5 * Milli(Metre));
    //ephemeris_->Prolong(JulianDate(2434378.5003725));  // Prolong to 1953-01-01
                                                       // 00:00:00 UTC. This
                                                       // *will* cause problems
                                                       // if not changed should
                                                       // more dates be added.
  }*/
  static void SetUpTestCase() {
    solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
    ephemeris_ = solar_system_1950_.MakeEphemeris(
        McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
        45 * Minute, 5 * Milli(Metre));
  }

  void CheckLunarUmbralEclipse(Instant const& current_time,
                               Sign const moon_offset_sign) {
    ephemeris_->Prolong(current_time);
    auto const sun = solar_system_1950_.massive_body(*ephemeris_, "Sun");
    auto const earth = solar_system_1950_.massive_body(*ephemeris_, "Earth");
    auto const moon = solar_system_1950_.massive_body(*ephemeris_, "Moon");

    auto const q_sun = ephemeris_->trajectory(sun)
                           ->EvaluatePosition(current_time, /*hint=*/nullptr);
    auto const q_moon = ephemeris_->trajectory(moon)
                            ->EvaluatePosition(current_time, /*hint=*/nullptr);
    auto const q_earth = ephemeris_->trajectory(earth)
                             ->EvaluatePosition(current_time, /*hint=*/nullptr);

    // MassiveBody will need radius information.  Or non-hardcoded data from
    // https://github.com/mockingbirdnest/Principia/blob/master/astronomy/gravity_model.proto.txt
    auto const r_sun = 696000.0 * Kilo(Metre);
    auto const r_earth = 6378.1363 * Kilo(Metre);
    auto const r_moon = 1738.0 * Kilo(Metre);

    // Checking body angles at the target time.
    // Angle formed by a right circular cone with sides defined by tangent lines
    // between Sun and Earth, and axis running through the centers of each.
    auto const umbral_half_aperture =
        ArcSin((r_sun - r_earth) / (q_sun - q_earth).Norm());
    auto const apex_of_moon_locus_at_umbral_contact =
        q_earth +
        Normalize(q_earth - q_sun) * (r_earth + moon_offset_sign * r_moon) /
            Sin(umbral_half_aperture);
    // Angle between Earth and Moon as seen at
    // apex_of_moon_locus_at_umbral_contact.
    auto const earth_moon_angle =
        ArcCos(InnerProduct(apex_of_moon_locus_at_umbral_contact - q_earth,
                            apex_of_moon_locus_at_umbral_contact - q_moon) /
               ((apex_of_moon_locus_at_umbral_contact - q_moon).Norm() *
                (apex_of_moon_locus_at_umbral_contact - q_earth).Norm()));
    // We are at the desired contact if the angle between Earth and Moon from
    // the apex of locus of the moon at that contact is the same value as the
    // half-aperture of the umbra (Earth-Sun cone).
    EXPECT_THAT(AbsoluteError(umbral_half_aperture, earth_moon_angle),
                AllOf(Lt(1.0 * Milli(Radian)), Gt(1.0 * Nano(Radian))))
        << NAMED(umbral_half_aperture) << ", " << NAMED(earth_moon_angle)
        << ", " << NAMED(current_time);
  }

  void CheckLunarPenumbralEclipse(Instant const& current_time,
                                  Sign const moon_offset_sign) {
    ephemeris_->Prolong(current_time);
    auto const sun = solar_system_1950_.massive_body(*ephemeris_, "Sun");
    auto const earth = solar_system_1950_.massive_body(*ephemeris_, "Earth");
    auto const moon = solar_system_1950_.massive_body(*ephemeris_, "Moon");

    auto const q_sun = ephemeris_->trajectory(sun)
                           ->EvaluatePosition(current_time, /*hint=*/nullptr);
    auto const q_moon = ephemeris_->trajectory(moon)
                            ->EvaluatePosition(current_time, /*hint=*/nullptr);
    auto const q_earth = ephemeris_->trajectory(earth)
                             ->EvaluatePosition(current_time, /*hint=*/nullptr);

    // MassiveBody will need radius information.  Or non-hardcoded data from
    // https://github.com/mockingbirdnest/Principia/blob/master/astronomy/gravity_model.proto.txt
    auto const r_sun = 696000.0 * Kilo(Metre);
    auto const r_earth = 6378.1363 * Kilo(Metre);
    auto const r_moon = 1738.0 * Kilo(Metre);

    auto const penumbral_half_aperture =
        ArcSin((r_sun + r_earth) / (q_sun - q_earth).Norm());
    auto const apex_of_moon_locus_at_penumbral_contact =
        q_earth +
        Normalize(q_sun - q_earth) * (r_earth + moon_offset_sign * r_moon) /
            Sin(penumbral_half_aperture);
    // Angle between Earth and Moon as seen at
    // apex_of_moon_locus_at_penumbral_contact.
    auto const earth_moon_angle =
        ArcCos(InnerProduct(apex_of_moon_locus_at_penumbral_contact - q_earth,
                            apex_of_moon_locus_at_penumbral_contact - q_moon) /
               ((apex_of_moon_locus_at_penumbral_contact - q_moon).Norm() *
                (apex_of_moon_locus_at_penumbral_contact - q_earth).Norm()));
    // We are at the desired contact if the angle between Earth and Moon from
    // the apex of locus of the moon at that contact is the same value as the
    // half-aperture of the penumbra.
    EXPECT_THAT(AbsoluteError(penumbral_half_aperture, earth_moon_angle),
                AllOf(Lt(1.0 * Milli(Radian)), Gt(1.0 * Nano(Radian))))
        << NAMED(penumbral_half_aperture) << ", " << NAMED(earth_moon_angle)
        << ", " << NAMED(current_time);
  }

  static SolarSystem<ICRFJ2000Equator> solar_system_1950_;
  static std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
};

#if !defined(_DEBUG)
TEST_F(EclipseTest, Year1950) {
  // Dates are TDB Julian Day for 1950-04-02.
  auto P1 = JulianDate(2433374.25788409);  // 18:10:49 UT
  auto U1 = JulianDate(2433374.29850909);  // 19:09:19
  auto U2 = JulianDate(2433374.354979);    // 20:30:38
  auto U3 = JulianDate(2433374.37367113);  // 20:57:33
  auto U4 = JulianDate(2433374.43016419);  // 22:18:54
  auto P4 = JulianDate(2433374.47075446);  // 23:17:21

  CheckLunarPenumbralEclipse(P1, U14);
  CheckLunarUmbralEclipse(U1, U14);
  CheckLunarUmbralEclipse(U2, U23);
  CheckLunarUmbralEclipse(U3, U23);
  CheckLunarUmbralEclipse(U4, U14);
  CheckLunarPenumbralEclipse(P4, U14);

  // Dates are TDB Julian Day for 1950-09-26.
  P1 = JulianDate(2433550.55712016);  // 01:21:43 UT
  U1 = JulianDate(2433550.60578913);  // 02:31:48
  U2 = JulianDate(2433550.66325441);  // 03:54:33
  U3 = JulianDate(2433550.69399515);  // 04:38:49
  U4 = JulianDate(2433550.75144885);  // 06:01:33
  P4 = JulianDate(2433550.800222);    // 07:11:47

  CheckLunarPenumbralEclipse(P1, U14);
  CheckLunarUmbralEclipse(U1, U14);
  CheckLunarUmbralEclipse(U2, U23);
  CheckLunarUmbralEclipse(U3, U23);
  CheckLunarUmbralEclipse(U4, U14);
  CheckLunarPenumbralEclipse(P4, U14);
}

TEST_F(EclipseTest, Year1951) {
  // Dates are TDB Julian Day for 1951-03-23.
  auto P1 = JulianDate(2433728.86842806);  // 08:50:50
  auto P4 = JulianDate(2433729.01725909);  // 12:24:19

  CheckLunarPenumbralEclipse(P1, U14);
  CheckLunarPenumbralEclipse(P4, U14);

  // Dates are TDB Julian Day for 1951-09-15.
  P1 = JulianDate(2433904.93736321);  // 10:29:16
  P4 = JulianDate(2433905.1002799);   // 14:23:52

  CheckLunarPenumbralEclipse(P1, U14);
  CheckLunarPenumbralEclipse(P4, U14);
}

TEST_F(EclipseTest, Year1952) {
  // Dates are TDB Julian Day for 1952-02-11 (or 10 for P1).
  auto P1 = JulianDate(2434053.42282623);  // P1 = 22:08:20 UT
  auto U1 = JulianDate(2434053.50334705);  // U1 = 00:04:17
  auto U4 = JulianDate(2434053.55203917);  // U4 = 01:14:24
  auto P4 = JulianDate(2434053.63249055);  // P4 = 03:10:15

  CheckLunarPenumbralEclipse(P1, U14);
  CheckLunarUmbralEclipse(U1, U14);
  CheckLunarUmbralEclipse(U4, U14);
  CheckLunarPenumbralEclipse(P4, U14);

  // Dates are TDB Julian Day for 1952-08-05.
  P1 = JulianDate(2434230.22830075);  // P1 = 17:28:13 UT
  U1 = JulianDate(2434230.27385631);  // U1 = 18:33:49
  U4 = JulianDate(2434230.37606695);  // U4 = 21:01:00
  P4 = JulianDate(2434230.42161093);  // P4 = 22:06:35

  CheckLunarPenumbralEclipse(P1, U14);
  CheckLunarUmbralEclipse(U1, U14);
  CheckLunarUmbralEclipse(U4, U14);
  CheckLunarPenumbralEclipse(P4, U14);

  // Later on for additional accuracy: 2 * ArcTan((x_norm_y -
  // y_normx).Norm(),(x_norm_y + y_norm_x).Norm())
  // x_norm_y = x * y.Norm() and y_norm_x = y * x.Norm()

  // Future is not now: check 2048-01-01 Lunar eclipse.
  /*P1 = JulianDate(2469076.66235167);  // 03:52:39 UT
  U1 = JulianDate(2469076.71279148);  // 05:05:17
  U2 = JulianDate(2469076.76776833);  // 06:24:27
  U3 = JulianDate(2469076.80661092);  // 07:20:23
  U4 = JulianDate(2469076.86158778);  // 08:39:33
  P4 = JulianDate(2469076.91195815);  // 09:52:05

  CheckLunarPenumbralEclipse(P1, U14);
  CheckLunarUmbralEclipse(U1, U14);
  CheckLunarUmbralEclipse(U2, U23);
  CheckLunarUmbralEclipse(U3, U23);
  CheckLunarUmbralEclipse(U4, U14);
  CheckLunarPenumbralEclipse(P4, U14); */
}
#endif

}  // namespace physics
}  // namespace principia
