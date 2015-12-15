#include "geometry/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/solar_system.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::JulianDate;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::ArcCos;
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

class EclipseTest : public testing::Test {
 protected:
  EclipseTest() {
    solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "initial_state_jd_2433282_500000000.proto.txt");
  }

  SolarSystem<ICRFJ2000Equator> solar_system_1950_;
};

TEST_F(EclipseTest, Dummy) {
  auto ephemeris = solar_system_1950_.MakeEphemeris(
      McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
      45 * Minute, 5 * Milli(Metre));

  ephemeris->Prolong(JulianDate(2433374.5));  // Prolong just past date of eclipse
                                              // (Eclipse was 1950-04-02 but JD
                                              // is 1950-04-03:00:00:00).
  // Pass body to Ephemeris.trajectory.
  auto const sun = solar_system_1950_.massive_body(*ephemeris, "Sun");
  auto const earth = solar_system_1950_.massive_body(*ephemeris, "Earth");
  auto const moon = solar_system_1950_.massive_body(*ephemeris, "Moon");

  // MassiveBody eventually needs radius information. Or non-hardcoded data from
  // https://github.com/mockingbirdnest/Principia/blob/master/astronomy/gravity_model.proto.txt
  auto const r_sun = 696000.0 * Kilo(Metre);
  auto const r_earth = 6378.1363 * Kilo(Metre);
  auto const r_moon = 1738.0 * Kilo(Metre);

  // Dates are TDB Julian Day for 1948-04-02.
  auto const P1 = JulianDate(2433374.25788409);  // 18:10:49 UT
  auto const U1 = JulianDate(2433374.29850909);  // 19:09:19
  auto const U2 = JulianDate(2433374.354979);    // 20:30:38
  auto const U3 = JulianDate(2433374.37367113);  // 20:57:33
  auto const U4 = JulianDate(2433374.43016419);  // 22:18:54
  auto const P4 = JulianDate(2433374.47075446);  // 23:17:21

  // Lunar eclipse -- U1 (time is start of umbral phase).
  auto q_sun = ephemeris->trajectory(sun)->EvaluatePosition(U1, nullptr);
  auto q_moon = ephemeris->trajectory(moon)->EvaluatePosition(U1, nullptr);
  auto q_earth = ephemeris->trajectory(earth)->EvaluatePosition(U1, nullptr);

  // Check body angles at target time.
  // Earth/Sun lineup.
  auto half_sun_earth_aperture = ArcSin((r_sun - r_earth) / (q_sun - q_earth).Norm());
  auto q_U14 =
      q_earth + Normalize(q_earth - q_sun) * (r_earth + r_moon) / Sin(half_sun_earth_aperture);
  // Earth/Moon lineup for U14.
  auto half_earth_moon_aperture = ArcCos(InnerProduct(q_U14 - q_earth, q_U14 - q_moon) /
                      ((q_U14 - q_moon).Norm() * (q_U14 - q_earth).Norm()));
  // Do the angles work?
  EXPECT_THAT(AbsoluteError(half_sun_earth_aperture, half_earth_moon_aperture),
              AllOf(Lt(1.0 * Milli(Radian)), Gt(1.0 * Nano(Radian)))) <<
              NAMED(half_sun_earth_aperture) << ", " << NAMED(half_earth_moon_aperture);

  // Lunar eclipse -- U2 (time is start of full eclipse).
  q_sun = ephemeris->trajectory(sun)->EvaluatePosition(U2, nullptr);
  q_moon = ephemeris->trajectory(moon)->EvaluatePosition(U2, nullptr);
  q_earth = ephemeris->trajectory(earth)->EvaluatePosition(U2, nullptr);

  // Check body angles at target time.
  // Earth/Sun lineup.
  half_sun_earth_aperture = ArcSin((r_sun - r_earth) / (q_sun - q_earth).Norm());
  auto q_U23 =
      q_earth + Normalize(q_earth - q_sun) * (r_earth - r_moon) / Sin(half_sun_earth_aperture);
  // Earth/Moon lineup for U23.
  half_earth_moon_aperture = ArcCos(InnerProduct(q_U23 - q_earth, q_U23 - q_moon) /
                     ((q_U23 - q_moon).Norm() * (q_U23 - q_earth).Norm()));
  // Do the angles work?
  EXPECT_THAT(AbsoluteError(half_sun_earth_aperture, half_earth_moon_aperture),
              AllOf(Lt(1.0 * Milli(Radian)), Gt(1.0 * Nano(Radian)))) <<
              NAMED(half_sun_earth_aperture) << ", " << NAMED(half_earth_moon_aperture);

  // Lunar eclipse -- U3 (time is end of full eclipse).
  q_sun = ephemeris->trajectory(sun)->EvaluatePosition(U3, nullptr);
  q_moon = ephemeris->trajectory(moon)->EvaluatePosition(U3, nullptr);
  q_earth = ephemeris->trajectory(earth)->EvaluatePosition(U3, nullptr);

  // Check body angles at target time.
  // Earth/Sun lineup.
  half_sun_earth_aperture = ArcSin((r_sun - r_earth) / (q_sun - q_earth).Norm());
  q_U23 =
      q_earth + Normalize(q_earth - q_sun) * (r_earth - r_moon) / Sin(half_sun_earth_aperture);
  // Earth/Moon lineup for U23.
  half_earth_moon_aperture = ArcCos(InnerProduct(q_U23 - q_earth, q_U23 - q_moon) /
                ((q_U23 - q_moon).Norm() * (q_U23 - q_earth).Norm()));
  // Do the angles work?
  EXPECT_THAT(AbsoluteError(half_sun_earth_aperture, half_earth_moon_aperture),
              AllOf(Lt(1.0 * Milli(Radian)), Gt(1.0 * Nano(Radian)))) <<
              NAMED(half_sun_earth_aperture) << ", " << NAMED(half_earth_moon_aperture);

  // Lunar eclipse -- U4 (time is end of umbral phase).
  q_sun = ephemeris->trajectory(sun)->EvaluatePosition(U4, nullptr);
  q_moon = ephemeris->trajectory(moon)->EvaluatePosition(U4, nullptr);
  q_earth = ephemeris->trajectory(earth)->EvaluatePosition(U4, nullptr);

  // Check body angles at target time.
  // Earth/Sun lineup.
  half_sun_earth_aperture = ArcSin((r_sun - r_earth) / (q_sun - q_earth).Norm());
  q_U14 =
      q_earth + Normalize(q_earth - q_sun) * (r_earth + r_moon) / Sin(half_sun_earth_aperture);
  // Earth/Moon lineup for U14.
  half_earth_moon_aperture = ArcCos(InnerProduct(q_U14 - q_earth, q_U14 - q_moon) /
                 ((q_U14 - q_moon).Norm() * (q_U14 - q_earth).Norm()));
  // Do the angles work?
  EXPECT_THAT(AbsoluteError(half_sun_earth_aperture, half_earth_moon_aperture),
              AllOf(Lt(1.0 * Milli(Radian)), Gt(1.0 * Nano(Radian)))) <<
              NAMED(half_sun_earth_aperture) << ", " << NAMED(half_earth_moon_aperture);

  // Later on for additional accuracy: 2 * ArcTan((x_norm_y -
  // y_normx).Norm(),(x_norm_y + y_norm_x).Norm())
  // x_norm_y = x * y.Norm() and y_norm_x = y * x.Norm()

  // Future: check 2048-01-01 Lunar eclipse.
  // P1 = 03:52:39 UT
  // U1 = 05:05:17 UT
  // U2 = 06:24:27
  // U3 = 07:20:23
  // U4 = 08:39:33
  // P4 = 09:52:05
  // etimes = {2469076.66235167, 2469076.71279148, 2469076.76776833,
  // 2469076.80661092, 2469076.86158778, 2469076.91195815};
}

}  // namespace physics
}  // namespace principia
