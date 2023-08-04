#include "ksp_plugin/flight_plan_optimizer.hpp"

#include "astronomy/date_time.hpp"
#include "astronomy/time_scales.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin_test/plugin_io.hpp"
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

TEST_F(FlightPlanOptimizerTest, Reach) {
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
  LOG(ERROR)<<flyby_time<<" "<<flyby_distance;
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:02:40"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(58591.4_(1) * Kilo(Metre)));

  FlightPlanOptimizer optimizer(&flight_plan);

  //auto const manœuvre0 = flight_plan.GetManœuvre(0);
  //CHECK_OK(optimizer.Optimize(/*index=*/0, moon, 1 * Metre / Second));

  //ComputeFlyby(flight_plan, moon, flyby_time, flyby_distance);
  //LOG(ERROR)<<flyby_time<<" "<<flyby_distance;
  //EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-26T21:00:36"_DateTime));
  //EXPECT_THAT(flyby_distance, IsNear(44949.7_(1) * Kilo(Metre)));

  //CHECK_OK(flight_plan.Replace(manœuvre0.burn(), /*index=*/0));

  //auto const manœuvre1 = flight_plan.GetManœuvre(1);
  //CHECK_OK(optimizer.Optimize(/*index=*/1, moon, 1 * Metre / Second));

  //ComputeFlyby(flight_plan, moon, flyby_time, flyby_distance);
  //LOG(ERROR)<<flyby_time<<" "<<flyby_distance;
  //EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:07:23"_DateTime));
  //EXPECT_THAT(flyby_distance, IsNear(57785.6_(1) * Kilo(Metre)));

  //CHECK_OK(flight_plan.Replace(manœuvre1.burn(), /*index=*/1));

  //auto const manœuvre2 = flight_plan.GetManœuvre(2);
  //CHECK_OK(optimizer.Optimize(/*index=*/2, moon, 1 * Metre / Second));

  //ComputeFlyby(flight_plan, moon, flyby_time, flyby_distance);
  //LOG(ERROR)<<flyby_time<<" "<<flyby_distance;
  //EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:02:26"_DateTime));
  //EXPECT_THAT(flyby_distance, IsNear(58664.1_(1) * Kilo(Metre)));

  //CHECK_OK(flight_plan.Replace(manœuvre2.burn(), /*index=*/2));

  //auto const manœuvre3 = flight_plan.GetManœuvre(3);
  //CHECK_OK(optimizer.Optimize(/*index=*/3, moon, 1 * Metre / Second));

  //ComputeFlyby(flight_plan, moon, flyby_time, flyby_distance);
  //LOG(ERROR)<<flyby_time<<" "<<flyby_distance;
  //EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-26T21:04:31"_DateTime));
  //EXPECT_THAT(flyby_distance, IsNear(48062.2_(1) * Kilo(Metre)));

  //CHECK_OK(flight_plan.Replace(manœuvre3.burn(), /*index=*/3));

  auto const manœuvre6 = flight_plan.GetManœuvre(6);
  CHECK_OK(optimizer.Optimize(/*index=*/6, moon, 1 * Metre / Second));

  EXPECT_EQ(0, flight_plan.number_of_anomalous_manœuvres());
  EXPECT_EQ(manœuvre6.initial_time(),
            flight_plan.GetManœuvre(6).initial_time());
  EXPECT_THAT((manœuvre6.Δv() - flight_plan.GetManœuvre(6).Δv()).Norm(),
              IsNear(1.0_(1) * Metre / Second));

  ComputeFlyby(flight_plan, moon, flyby_time, flyby_distance);
  LOG(ERROR)<<flyby_time<<" "<<flyby_distance;
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:22:55"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(24788.5_(1) * Kilo(Metre)));

  CHECK_OK(flight_plan.Replace(manœuvre6.burn(), /*index=*/6));

  auto const manœuvre7 = flight_plan.GetManœuvre(7);
  EXPECT_THAT(optimizer.Optimize(/*index=*/7, moon, 1 * Metre / Second),
              StatusIs(termination_condition::VanishingStepSize));

  // We cannot compute flybys because the flight plan basically goes through the
  // centre of the Moon.
  EXPECT_EQ(8, flight_plan.number_of_anomalous_manœuvres());
  EXPECT_THAT(
      manœuvre7.initial_time() - flight_plan.GetManœuvre(7).initial_time(),
      IsNear(-3.4_(1) * Milli(Second)));
  EXPECT_THAT(
      (manœuvre7.Δv() - flight_plan.GetManœuvre(7).Δv()).Norm(),
      IsNear(61.9_(1) * Metre / Second));

  CHECK_OK(flight_plan.Replace(manœuvre7.burn(), /*index=*/7));
}

}  // namespace ksp_plugin
}  // namespace principia
