#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "physics/solar_system.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;

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
    
    TEST_F(EclipseTest, Dummy){
        // ...
    };
        
} // physics

} // principia
