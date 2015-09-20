#include "physics/solar_system_factory.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace physics {

class SolarSystemFactoryTest : public ::testing::Test {
 protected:
  SolarSystemFactory factory_;
};

TEST_F(SolarSystemFactoryTest, Parse) {
  factory_.Initialize(SOLUTION_DIR "astronomy\\gravity_model.proto.txt",
                      SOLUTION_DIR "astronomy\\initial_state.proto.txt");
}

}  // namespace physics
}  // namespace principia
