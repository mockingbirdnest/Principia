#include "ksp_plugin/flight_plan.hpp"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {

namespace ksp_plugin {

class FlightPlanTest : public testing::Test {

  FlightPlanTest() {
    std::vector<not_null<std::unique_ptr<MassiveBody>>> bodies;
    bodies.emplace_back(1 * Pow<3>(Metre) / Pow<2>(Second));
    DegreesOfFreedom<Barycentric> const dof = {Barycentric::origin,
                                               Velocity<Barycentric>()};
    Instant const t0;
    ephemeris_ = Ephemeris<Barycentric>(
        bodies,
        /*initial_state=*/{dof},
        /*initial_time=*/t0,
        integrators::McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
        /*step=*/1 * Second,
        /*fitting_tolerance=*/1 * Milli(Metre));
    root_.Append(t0, )
  }

  std::unique_ptr<Ephemeris<Barycentric>> ephemeris_;
  DiscreteTrajectory<Barycentric> root_;
  std::unique_ptr<FlightPlan> flight_plan_;
}

TEST_F(FlightPlanTest, Append) {
  
}

}  // namespace ksp_plugin
}  // namespace principia
