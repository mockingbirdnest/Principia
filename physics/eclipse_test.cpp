#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "geometry/epoch.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::JulianDate;
using integrators::McLachlanAtela1992Order5Optimal;

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
        auto ephemeris = solar_system_1950_.MakeEphemeris(McLachlanAtela1992Order5Optimal<Position<ICRFJ000Equator>>(),
		45 * Minute,
		5 * Milli(Metre));

		ephemeris->Prolong(JulianDate(2433374.5));
		// Prolong until date of eclipse
		// pass body to Ephemeris.trajectory
    };
        
} // physics

} // principia
