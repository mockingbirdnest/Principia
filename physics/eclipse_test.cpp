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

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::JulianDate;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::si::Minute;
using quantities::si::Metre;
using quantities::si::Milli;

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

    // Get positions/trajectories/ephemeres for bodies
    // Dates will have to be for U1, etc. (and all will have to be generated)
    auto some_instant = JulianDate(2433374.5);
    auto q_sun = ephemeris->trajectory(sun)->EvaluatePosition(some_instant, nullptr);
    ephemeris->trajectory(earth)->EvaluatePosition(some_instant, nullptr);
    ephemeris->trajectory(moon)->EvaluatePosition(some_instant, nullptr);

    // check body angles at target times
    // ArcTan(r_sun/Norm(q_moon - q_sun));
    // ArcTan(r_earth/Norm(q_moon - q_earth));
    LOG(ERROR) << ArcTan(1.0);
    // Future: check 2048-01-01 Lunar eclipse
  };
        
} // physics

} // principia
