
#include "ksp_plugin/interface.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin_test/mock_planetarium.hpp"
#include "ksp_plugin_test/mock_plugin.hpp"
#include "ksp_plugin_test/mock_renderer.hpp"
#include "ksp_plugin_test/mock_vessel.hpp"

namespace principia {
namespace interface {

using base::make_not_null_unique;
using ksp_plugin::FlightPlan;
using ksp_plugin::MockPlanetarium;
using ksp_plugin::MockPlugin;
using ksp_plugin::MockVessel;
using ::testing::StrictMock;

class InterfaceExternalTest : public ::testing::Test {
 protected:
  InterfaceExternalTest() : flight_plan_(1 * Tonne, t0_,) {
    ON_CALL()
  }

  MockPlugin const plugin_;
  MockVessel const vessel_;
  Instant const t0_;
  std::unique_ptr<Ephemeris
  std::unique_ptr<FlightPlan> flight_plan_;
};

TEST_F(InterfaceExternalTest, GetNearestPlannedCoastDegreesOfFreedom) {
  
}

}  // namespace interface
}  // namespace principia
