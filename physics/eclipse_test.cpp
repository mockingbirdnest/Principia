#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "geometry/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/elementary_functions.hpp"
#include "geometry/grassmann.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::JulianDate;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::si::Minute;
using quantities::si::Metre;
using quantities::si::Nano;
using quantities::si::Milli;
using quantities::si::Kilo;
using quantities::ArcCos; // This feels very cargocultish, especially since I don't need it for ArcSin.
using testing_utilities::AbsoluteError;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::AllOf;

namespace physics {

  class EclipseTest : public testing::Test {
  protected:
    EclipseTest() {
      solar_system_1950_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model.proto.txt", 
        SOLUTION_DIR / "astronomy" / "initial_state_jd_2433282_500000000.proto.txt");
    }

    SolarSystem<ICRFJ2000Equator> solar_system_1950_;
  };
    
  TEST_F(EclipseTest, Dummy) {
    auto ephemeris = solar_system_1950_.MakeEphemeris(McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
    45 * Minute,
    5 * Milli(Metre));

    ephemeris->Prolong(JulianDate(2433374.5)); // Prolong until date of eclipse (Eclipse was 1950-04-02 but JD is 1950-04-03:00:00:00)
    // pass body to Ephemeris.trajectory
    auto sun = ephemeris->bodies()[solar_system_1950_.index("Sun")];
    auto earth = ephemeris->bodies()[solar_system_1950_.index("Earth")];
    auto moon = ephemeris->bodies()[solar_system_1950_.index("Moon")];
    
    // Massive_body eventually needs radius information. Or non-hardcoded data pulled from https://github.com/mockingbirdnest/Principia/blob/master/astronomy/gravity_model.proto.txt
    auto r_sun = 696000.0 * Kilo(Metre);
    auto r_earth = 6378.1363 * Kilo(Metre);
    auto r_moon = 1738.0 * Kilo(Metre);

    // Get positions/trajectories/ephemeres for bodies
    // Dates will have to be TDB Julian Day for U1, etc. (and all will have to be generated)
    // P1 = 18:10:49 UT
    // U1 = 19:09:19
    // U2 = 20:30:38
    // U3 = 20:57:33
    // U4 = 22:18:54
    // P4 = 23:17:21
    // 2433374.25788409, 2433374.29850909, 2433374.354979, 2433374.37367113, 2433374.43016419, 2433374.47075446
    double etimes[6] = {2433374.25788409, 2433374.29850909, 2433374.354979, 2433374.37367113, 2433374.43016419, 2433374.47075446};
    for(int i = 0; i <= 5; i++) {
      auto some_instant = JulianDate(etimes[i]);
      auto q_sun = ephemeris->trajectory(sun)->EvaluatePosition(some_instant, nullptr);
      auto q_moon = ephemeris->trajectory(moon)->EvaluatePosition(some_instant, nullptr);
      auto q_earth = ephemeris->trajectory(earth)->EvaluatePosition(some_instant, nullptr);

      // check body angles at target times
      // Lunar eclipse
      // Earth/Sun lineup
      auto alpha = ArcSin((r_sun - r_earth)/(q_sun - q_earth).Norm());
      auto q_U23 = q_earth + Normalize(q_earth - q_sun) * (r_earth - r_moon) / Sin(alpha);
      auto q_U14 = q_earth + Normalize(q_earth - q_sun) * (r_earth + r_moon) / Sin(alpha);
      // Earth/Moon lineup
      auto beta = ArcCos(InnerProduct(q_U23 - q_earth, q_U23 - q_moon) / ((q_U23 - q_moon).Norm() * (q_U23 - q_earth).Norm())); // Angle for U23
      auto gamma = ArcCos(InnerProduct(q_U14 - q_earth, q_U14 - q_moon) / ((q_U14 - q_moon).Norm() * (q_U14 - q_earth).Norm())); // Angle for U14
      // Still need to compare angles
      LOG(ERROR) << alpha;
      LOG(ERROR) << beta;
      LOG(ERROR) << gamma;

      EXPECT_THAT(AbsoluteError(alpha, beta), AllOf(Lt(1 * Milli(Radian)), Gt(1 * Nano(Radian)))); // Should work for 3rd & 4th times
      EXPECT_THAT(AbsoluteError(alpha, gamma), AllOf(Lt(1 * Milli(Radian)), Gt(1 * Nano(Radian)))); // Should work for 2nd & 5th times
      }
    
    // Later on for additional accuracy: 2 * ArcTan((x_norm_y - y_normx).Norm(),(x_norm_y + y_norm_x).Norm())
    // x_norm_y = x * y.Norm() and y_norm_x = y * x.Norm()

    // Future: check 2048-01-01 Lunar eclipse
    // P1 = 03:52:39 UT
    // U1 = 05:05:17 UT
    // U2 = 06:24:27
    // U3 = 07:20:23
    // U4 = 08:39:33
    // P4 = 09:52:05
    // 2469076.66235167  2469076.71279148  2469076.76776833  2469076.80661092  2469076.86158778  2469076.91195815
    // etimes = {2469076.66235167, 2469076.71279148, 2469076.76776833, 2469076.80661092, 2469076.86158778, 2469076.91195815};
  };
        
} // physics

} // principia
