#include "ksp_plugin/flight_plan_optimizer.hpp"

#include <utility>
#include <vector>

#include "astronomy/date_time.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin_test/plugin_io.hpp"
#include "physics/apsides.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::Eq;
using ::testing::ResultOf;
using namespace principia::astronomy::_date_time;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_flight_plan_optimizer;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin_test::_plugin_io;
using namespace principia::physics::_apsides;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;

class FlightPlanOptimizerTest : public testing::Test {
 protected:
  FlightPlanOptimizerTest() {
    google::SetStderrLogging(google::INFO);
  }

  ~FlightPlanOptimizerTest() override {
    google::SetStderrLogging(FLAGS_stderrthreshold);
  }

  static Celestial const& FindCelestialByName(std::string_view const name,
                                              Plugin const& plugin) {
    for (Index index = 0; plugin.HasCelestial(index); ++index) {
      Celestial const& celestial = plugin.GetCelestial(index);
      if (celestial.body()->name() == name) {
        return celestial;
      }
    }
    LOG(FATAL) << "No celestial named " << name;
  }

  static void ComputeFlyby(FlightPlan const& flight_plan,
                           Celestial const& celestial,
                           Instant& flyby_time,
                           Length& flyby_distance) {
    auto const& celestial_trajectory = celestial.trajectory();
    auto const& flight_plan_trajectory = flight_plan.GetAllSegments();
    DiscreteTrajectory<Barycentric> apoapsides;
    DiscreteTrajectory<Barycentric> periapsides;
    ComputeApsides(celestial_trajectory,
                   flight_plan_trajectory,
                   flight_plan_trajectory.begin(),
                   flight_plan_trajectory.end(),
                   /*max_points=*/100,
                   apoapsides,
                   periapsides);
    auto const radius = celestial.body()->mean_radius();
    for (const auto [time, _] : periapsides) {
      Length const periapsis_distance =
          (celestial_trajectory.EvaluatePosition(time) -
           flight_plan_trajectory.EvaluatePosition(time))
              .Norm();
      if (periapsis_distance < 50 * radius) {
        flyby_time = time;
        flyby_distance = periapsis_distance;
      }
    }
  }
};

// This test tweaks the burns of the flight plan that precede the Moon flyby, in
// order to collide head on with the Moon.  This test takes about 10 minutes.
TEST_F(FlightPlanOptimizerTest, DISABLED_ReachTheMoon) {
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      SOLUTION_DIR / "ksp_plugin_test" / "saves" / "3072.proto.b64",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");

  auto const ifnity = plugin->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");
  EXPECT_THAT(ifnity->name(), Eq("IFNITY-5.2"));
  EXPECT_THAT(TTSecond(ifnity->trajectory().front().time),
              Eq("1970-08-14T08:03:46"_DateTime));
  EXPECT_THAT(TTSecond(ifnity->psychohistory()->back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  ASSERT_TRUE(ifnity->has_flight_plan());
  ifnity->ReadFlightPlanFromMessage();
  FlightPlan& flight_plan = ifnity->flight_plan();
  EXPECT_THAT(
      flight_plan.adaptive_step_parameters().length_integration_tolerance(),
      Eq(1 * Metre));
  EXPECT_THAT(flight_plan.adaptive_step_parameters().max_steps(), Eq(16'000));
  EXPECT_THAT(flight_plan.number_of_manœuvres(), Eq(16));
  std::vector<std::pair<DateTime, Speed>> manœuvre_ignition_tt_seconds_and_Δvs;
  for (int i = 0; i < flight_plan.number_of_manœuvres(); ++i) {
    manœuvre_ignition_tt_seconds_and_Δvs.emplace_back(
        TTSecond(flight_plan.GetManœuvre(i).initial_time()),
        flight_plan.GetManœuvre(i).Δv().Norm());
  }

  Celestial const& moon = FindCelestialByName("Moon", *plugin);
  Instant flyby_time;
  Length flyby_distance;
  ComputeFlyby(flight_plan, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:02:40"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(58591.4_(1) * Kilo(Metre)));

  FlightPlanOptimizer optimizer(&flight_plan);

  // In the code below we cannot compute flybys because the flight plan
  // basically goes through the centre of the Moon.

  LOG(INFO) << "Optimizing manœuvre 5";
  auto const manœuvre5 = flight_plan.GetManœuvre(5);
  EXPECT_THAT(optimizer.Optimize(/*index=*/5, moon, 1 * Milli(Metre) / Second),
              StatusIs(termination_condition::VanishingStepSize));

  EXPECT_EQ(8, flight_plan.number_of_anomalous_manœuvres());
  EXPECT_THAT(
      manœuvre5.initial_time() - flight_plan.GetManœuvre(5).initial_time(),
      IsNear(-10.3_(1) * Micro(Second)));
  EXPECT_THAT(
      (manœuvre5.Δv() - flight_plan.GetManœuvre(5).Δv()).Norm(),
      IsNear(1.09_(1) * Metre / Second));

  CHECK_OK(flight_plan.Replace(manœuvre5.burn(), /*index=*/5));

  LOG(INFO) << "Optimizing manœuvre 6";
  auto const manœuvre6 = flight_plan.GetManœuvre(6);
  EXPECT_THAT(optimizer.Optimize(/*index=*/6, moon, 1 * Milli(Metre) / Second),
              StatusIs(termination_condition::VanishingStepSize));

  EXPECT_EQ(8, flight_plan.number_of_anomalous_manœuvres());
  EXPECT_EQ(manœuvre6.initial_time(),
            flight_plan.GetManœuvre(6).initial_time());
  EXPECT_THAT((manœuvre6.Δv() - flight_plan.GetManœuvre(6).Δv()).Norm(),
              IsNear(1.29_(1) * Metre / Second));

  CHECK_OK(flight_plan.Replace(manœuvre6.burn(), /*index=*/6));

  LOG(INFO) << "Optimizing manœuvre 7";
  auto const manœuvre7 = flight_plan.GetManœuvre(7);
  EXPECT_THAT(optimizer.Optimize(/*index=*/7, moon, 1 * Milli(Metre) / Second),
              StatusIs(termination_condition::VanishingStepSize));

  EXPECT_EQ(8, flight_plan.number_of_anomalous_manœuvres());
  EXPECT_THAT(
      manœuvre7.initial_time() - flight_plan.GetManœuvre(7).initial_time(),
      IsNear(-4.9_(1) * Milli(Second)));
  EXPECT_THAT(
      (manœuvre7.Δv() - flight_plan.GetManœuvre(7).Δv()).Norm(),
      IsNear(62.3_(1) * Metre / Second));

  CHECK_OK(flight_plan.Replace(manœuvre7.burn(), /*index=*/7));
}

}  // namespace ksp_plugin
}  // namespace principia
