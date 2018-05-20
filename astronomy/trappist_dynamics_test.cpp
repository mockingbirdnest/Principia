
#include "astronomy/frames.hpp"
#include "gtest/gtest.h"
#include "physics/solar_system.hpp"

namespace principia {

using physics::SolarSystem;

namespace astronomy {

class TrappistDynamicsTest : public ::testing::Test {
 protected:
};

TEST_F(TrappistDynamicsTest, Smoke) {
  SolarSystem<Trappist> trappist_system(
      SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "trappist_initial_state_jd_2457282_805700000.proto.txt");
}

}  // namespace astronomy
}  // namespace principia
