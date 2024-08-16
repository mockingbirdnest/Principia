#include "ksp_plugin/flight_plan_optimizer.hpp"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "astronomy/date_time.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin_test/plugin_io.hpp"
#include "numerics/transposed_view.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::AnyOf;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Le;
using ::testing::Matcher;
using ::testing::ResultOf;
using namespace principia::astronomy::_date_time;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_flight_plan_optimizer;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin_test::_plugin_io;
using namespace principia::numerics::_transposed_view;
using namespace principia::physics::_apsides;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_numerics_matchers;

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

  static void ComputeFlyby(
      FlightPlan const& flight_plan,
      Celestial const& celestial,
      Instant& flyby_time,
      DegreesOfFreedom<Barycentric>& flyby_degrees_of_freedom) {
    auto const& celestial_trajectory = celestial.trajectory();
    auto const& flight_plan_trajectory = flight_plan.GetAllSegments();
    DiscreteTrajectory<Barycentric> apoapsides;
    DiscreteTrajectory<Barycentric> periapsides;
    ComputeApsides(celestial_trajectory,
                   flight_plan_trajectory,
                   flight_plan_trajectory.begin(),
                   flight_plan_trajectory.end(),
                   /*t_max=*/InfiniteFuture,
                   /*max_points=*/100,
                   apoapsides,
                   periapsides);
    auto const radius = celestial.body()->mean_radius();
    for (auto const& [time, degrees_of_freedom] : periapsides) {
      Length const periapsis_distance =
          (celestial_trajectory.EvaluatePosition(time) -
           degrees_of_freedom.position())
              .Norm();
      if (periapsis_distance < 50 * radius) {
        flyby_time = time;
        flyby_degrees_of_freedom = degrees_of_freedom;
      }
    }
  }

  static void ComputeFlyby(FlightPlan const& flight_plan,
                           Celestial const& celestial,
                           Instant& flyby_time,
                           Length& flyby_distance) {
    auto const& celestial_trajectory = celestial.trajectory();
    DegreesOfFreedom<Barycentric> flyby_degrees_of_freedom(
        Barycentric::origin, Barycentric::unmoving);
    ComputeFlyby(flight_plan, celestial, flyby_time, flyby_degrees_of_freedom);
    flyby_distance = (celestial_trajectory.EvaluatePosition(flyby_time) -
                      flyby_degrees_of_freedom.position())
                         .Norm();
  }

  static void ComputeFlyby(FlightPlan const& flight_plan,
                           Celestial const& celestial,
                           NavigationFrame const& frame,
                           Instant& flyby_time,
                           Angle& flyby_inclination) {
    DegreesOfFreedom<Barycentric> flyby_degrees_of_freedom(
        Barycentric::origin, Barycentric::unmoving);
    ComputeFlyby(flight_plan, celestial, flyby_time, flyby_degrees_of_freedom);
    auto const navigation_degrees_of_freedom =
        frame.ToThisFrameAtTime(flyby_time)(flyby_degrees_of_freedom);
    auto const r =
        navigation_degrees_of_freedom.position() - Navigation::origin;
    auto const& v = navigation_degrees_of_freedom.velocity();
    flyby_inclination =
        AngleBetween(Wedge(r, v), Bivector<double, Navigation>({0, 0, 1}));
  }

  void ReadReachFlightPlan() {
    plugin_ = ReadPluginFromFile(
        SOLUTION_DIR / "ksp_plugin_test" / "saves" / "3072.proto.b64",
        /*compressor=*/"gipfeli",
        /*decoder=*/"base64");

    auto const ifnity =
        plugin_->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");
    EXPECT_THAT(ifnity->name(), Eq("IFNITY-5.2"));
    EXPECT_THAT(TTSecond(ifnity->trajectory().front().time),
                Eq("1970-08-14T08:03:46"_DateTime));
    EXPECT_THAT(TTSecond(ifnity->psychohistory()->back().time),
                Eq("1970-08-14T08:47:05"_DateTime));
    ASSERT_TRUE(ifnity->has_flight_plan());
    ifnity->ReadFlightPlanFromMessage();
    flight_plan_ = &ifnity->flight_plan();
    EXPECT_THAT(
        flight_plan_->adaptive_step_parameters().length_integration_tolerance(),
        Eq(1 * Metre));
    EXPECT_THAT(flight_plan_->adaptive_step_parameters().max_steps(),
                Eq(16'000));
    EXPECT_THAT(flight_plan_->number_of_manœuvres(), Eq(16));
    std::vector<std::pair<DateTime, Speed>>
        manœuvre_ignition_tt_seconds_and_Δvs;
    for (int i = 0; i < flight_plan_->number_of_manœuvres(); ++i) {
      manœuvre_ignition_tt_seconds_and_Δvs.emplace_back(
          TTSecond(flight_plan_->GetManœuvre(i).initial_time()),
          flight_plan_->GetManœuvre(i).Δv().Norm());
    }
  }

  std::unique_ptr<Plugin const> plugin_;
  FlightPlan* flight_plan_;
};

// This test tweaks the burns of the flight plan that precede the Moon flyby, in
// order to collide head on with the Moon.  This test takes about 10 minutes.
TEST_F(FlightPlanOptimizerTest, DISABLED_ReachTheMoon) {
  ReadReachFlightPlan();

  Celestial const& moon = FindCelestialByName("Moon", *plugin_);
  Instant flyby_time;
  Length flyby_distance;
  ComputeFlyby(*flight_plan_, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:02:40"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(58591.4_(1) * Kilo(Metre)));

  std::int64_t number_of_evaluations = 0;
  FlightPlanOptimizer optimizer(
      flight_plan_,
      FlightPlanOptimizer::ForCelestialCentre(&moon),
      [&number_of_evaluations](FlightPlan const&) { ++number_of_evaluations; });

  // In the code below we cannot compute flybys because the flight plan
  // basically goes through the centre of the Moon.

  LOG(INFO) << "Optimizing manœuvre 5";
  auto const manœuvre5 = flight_plan_->GetManœuvre(5);
  EXPECT_THAT(optimizer.Optimize(/*index=*/5, 1 * Milli(Metre) / Second),
              StatusIs(termination_condition::VanishingStepSize));

  EXPECT_EQ(8, flight_plan_->number_of_anomalous_manœuvres());
  EXPECT_THAT(
      manœuvre5.initial_time() - flight_plan_->GetManœuvre(5).initial_time(),
      IsNear(7.39_(1) * Micro(Second)));
  EXPECT_THAT(
      (manœuvre5.Δv() - flight_plan_->GetManœuvre(5).Δv()).Norm(),
      IsNear(1.054_(1) * Metre / Second));
  EXPECT_EQ(113, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre5.burn(), /*index=*/5));

  LOG(INFO) << "Optimizing manœuvre 6";
  auto const manœuvre6 = flight_plan_->GetManœuvre(6);
  EXPECT_THAT(optimizer.Optimize(/*index=*/6, 1 * Milli(Metre) / Second),
              StatusIs(termination_condition::VanishingStepSize));

  EXPECT_EQ(8, flight_plan_->number_of_anomalous_manœuvres());
  EXPECT_THAT(
      manœuvre6.initial_time() - flight_plan_->GetManœuvre(6).initial_time(),
      IsNear(12.1_(1) * Micro(Second)));
  EXPECT_THAT((manœuvre6.Δv() - flight_plan_->GetManœuvre(6).Δv()).Norm(),
              IsNear(1.292_(1) * Metre / Second));
  EXPECT_EQ(124, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre6.burn(), /*index=*/6));

  LOG(INFO) << "Optimizing manœuvre 7";
  auto const manœuvre7 = flight_plan_->GetManœuvre(7);
  EXPECT_THAT(optimizer.Optimize(/*index=*/7, 1 * Milli(Metre) / Second),
              StatusIs(termination_condition::VanishingStepSize));

  EXPECT_EQ(8, flight_plan_->number_of_anomalous_manœuvres());
  EXPECT_THAT(
      manœuvre7.initial_time() - flight_plan_->GetManœuvre(7).initial_time(),
      IsNear(-77_(1) * Milli(Second)));
  EXPECT_THAT(
      (manœuvre7.Δv() - flight_plan_->GetManœuvre(7).Δv()).Norm(),
      IsNear(61.9_(1) * Metre / Second));
  EXPECT_EQ(102, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre7.burn(), /*index=*/7));
}

// This test tweaks the burns of the flight plan that precede the Moon flyby, in
// order to fly at (approximately) 2000 km from the centre of the Moon.  This
// test takes about 10 minutes.
TEST_F(FlightPlanOptimizerTest, DISABLED_GrazeTheMoon) {
  ReadReachFlightPlan();

  Celestial const& moon = FindCelestialByName("Moon", *plugin_);
  Instant flyby_time;
  Length flyby_distance;
  ComputeFlyby(*flight_plan_, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:02:40"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(58591.4_(1) * Kilo(Metre)));

  std::int64_t number_of_evaluations = 0;
  FlightPlanOptimizer optimizer(
      flight_plan_,
      FlightPlanOptimizer::ForCelestialDistance(&moon, 2000 * Kilo(Metre)),
      [&number_of_evaluations](FlightPlan const&) { ++number_of_evaluations; });

  LOG(INFO) << "Optimizing manœuvre 5";
  auto const manœuvre5 = flight_plan_->GetManœuvre(5);
  EXPECT_OK(optimizer.Optimize(/*index=*/5, 1 * Milli(Metre) / Second));

  EXPECT_THAT(
      manœuvre5.initial_time() - flight_plan_->GetManœuvre(5).initial_time(),
      IsNear(43.8_(1) * Micro(Second)));
  EXPECT_THAT(
      (manœuvre5.Δv() - flight_plan_->GetManœuvre(5).Δv()).Norm(),
      IsNear(1.116_(1) * Metre / Second));

  ComputeFlyby(*flight_plan_, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:24:00"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(2255.3_(1) * Kilo(Metre)));
  EXPECT_EQ(79, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre5.burn(), /*index=*/5));

  LOG(INFO) << "Optimizing manœuvre 6";
  auto const manœuvre6 = flight_plan_->GetManœuvre(6);
  EXPECT_OK(optimizer.Optimize(/*index=*/6, 1 * Milli(Metre) / Second));

  EXPECT_THAT(
      manœuvre6.initial_time() - flight_plan_->GetManœuvre(6).initial_time(),
      IsNear(0.954_(1) * Micro(Second)));
  EXPECT_THAT((manœuvre6.Δv() - flight_plan_->GetManœuvre(6).Δv()).Norm(),
              IsNear(1.312_(1) * Metre / Second));

  ComputeFlyby(*flight_plan_, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:16:41"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(2001.4_(1) * Kilo(Metre)));
  EXPECT_EQ(72, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre6.burn(), /*index=*/6));

  LOG(INFO) << "Optimizing manœuvre 7";
  auto const manœuvre7 = flight_plan_->GetManœuvre(7);
  EXPECT_OK(optimizer.Optimize(/*index=*/7, 1 * Milli(Metre) / Second));

  EXPECT_THAT(
      manœuvre7.initial_time() - flight_plan_->GetManœuvre(7).initial_time(),
      IsNear(2.1_(1) * Milli(Second)));
  EXPECT_THAT(
      (manœuvre7.Δv() - flight_plan_->GetManœuvre(7).Δv()).Norm(),
      IsNear(59.3_(1) * Metre / Second));

  ComputeFlyby(*flight_plan_, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:14:59"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(2000.2_(1) * Kilo(Metre)));
  EXPECT_EQ(86, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre7.burn(), /*index=*/7));
}

TEST_F(FlightPlanOptimizerTest, DISABLED_GrindsToAHalt) {
  ReadReachFlightPlan();

  std::int64_t number_of_evaluations = 0;
  FlightPlanOptimizer optimizer(
      flight_plan_,
      FlightPlanOptimizer::ForΔv(),
      [&number_of_evaluations](FlightPlan const&) { ++number_of_evaluations; });

  LOG(INFO) << "Optimizing manœuvre 5";
  auto const manœuvre5 = flight_plan_->GetManœuvre(5);
  EXPECT_OK(optimizer.Optimize(/*index=*/5, 1 * Micro(Metre) / Second));

  // The initial time doesn't change.
  EXPECT_THAT(
      manœuvre5.initial_time() - flight_plan_->GetManœuvre(5).initial_time(),
      Eq(0 * Second));
  EXPECT_THAT(flight_plan_->GetManœuvre(5).Δv().Norm(),
              IsNear(4.3e-19_(1) * Metre / Second));
  EXPECT_EQ(0, number_of_evaluations);
}

TEST_F(FlightPlanOptimizerTest, DISABLED_PoleTheMoon) {
  ReadReachFlightPlan();

  Celestial const& moon = FindCelestialByName("Moon", *plugin_);
  auto const moon_index = plugin_->CelestialIndexOfBody(*moon.body());
  auto const moon_frame =
      plugin_->NewBodyCentredNonRotatingNavigationFrame(moon_index);
  Instant flyby_time;
  Angle flyby_inclination;
  ComputeFlyby(*flight_plan_, moon, *moon_frame, flyby_time, flyby_inclination);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:02:40"_DateTime));
  EXPECT_THAT(flyby_inclination, IsNear(76.32_(1) * Degree));

  std::int64_t number_of_evaluations = 0;
  FlightPlanOptimizer optimizer(
      flight_plan_,
      FlightPlanOptimizer::ForInclination(
          &moon,
          [this, moon_index]() {
            return plugin_->NewBodyCentredNonRotatingNavigationFrame(
                moon_index);
          },
          90 * Degree),
      [&number_of_evaluations](FlightPlan const&) { ++number_of_evaluations; });

  LOG(INFO) << "Optimizing manœuvre 5";
  auto const manœuvre5 = flight_plan_->GetManœuvre(5);
  EXPECT_OK(optimizer.Optimize(/*index=*/5, 1 * Milli(Metre) / Second));

  EXPECT_THAT(
      manœuvre5.initial_time() - flight_plan_->GetManœuvre(5).initial_time(),
      IsNear(0.715_(1) * Micro(Second)));
  EXPECT_THAT(
      (manœuvre5.Δv() - flight_plan_->GetManœuvre(5).Δv()).Norm(),
      IsNear(1.615_(1) * Centi(Metre) / Second));

  ComputeFlyby(*flight_plan_, moon, *moon_frame, flyby_time, flyby_inclination);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:09:20"_DateTime));
  EXPECT_THAT(flyby_inclination, IsNear(90.23_(1) * Degree));
  EXPECT_EQ(25, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre5.burn(), /*index=*/5));

  LOG(INFO) << "Optimizing manœuvre 6";
  auto const manœuvre6 = flight_plan_->GetManœuvre(6);
  EXPECT_OK(optimizer.Optimize(/*index=*/6, 1 * Milli(Metre) / Second));

  EXPECT_THAT(
      manœuvre6.initial_time() - flight_plan_->GetManœuvre(6).initial_time(),
      Eq(0 * Second));
  EXPECT_THAT((manœuvre6.Δv() - flight_plan_->GetManœuvre(6).Δv()).Norm(),
              IsNear(0.257_(1) * Metre / Second));

  ComputeFlyby(*flight_plan_, moon, *moon_frame, flyby_time, flyby_inclination);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:09:47"_DateTime));
  EXPECT_THAT(flyby_inclination, IsNear(89.98_(1) * Degree));
  EXPECT_EQ(34, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre6.burn(), /*index=*/6));

  LOG(INFO) << "Optimizing manœuvre 7";
  auto const manœuvre7 = flight_plan_->GetManœuvre(7);
  EXPECT_OK(optimizer.Optimize(/*index=*/7, 1 * Milli(Metre) / Second));

  EXPECT_THAT(
      manœuvre7.initial_time() - flight_plan_->GetManœuvre(7).initial_time(),
      IsNear(0.46_(1) * Milli(Second)));
  EXPECT_THAT(
      (manœuvre7.Δv() - flight_plan_->GetManœuvre(7).Δv()).Norm(),
      IsNear(12.7_(1) * Metre / Second));

  ComputeFlyby(*flight_plan_, moon, *moon_frame, flyby_time, flyby_inclination);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:08:39"_DateTime));
  EXPECT_THAT(flyby_inclination, IsNear(90.00_(1) * Degree));
  EXPECT_EQ(54, number_of_evaluations);
  number_of_evaluations = 0;

  CHECK_OK(flight_plan_->Replace(manœuvre7.burn(), /*index=*/7));
}

TEST_F(FlightPlanOptimizerTest, DISABLED_Combined) {
  ReadReachFlightPlan();

  Celestial const& moon = FindCelestialByName("Moon", *plugin_);
  Instant flyby_time;
  Length flyby_distance;
  ComputeFlyby(*flight_plan_, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:02:40"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(58591.4_(1) * Kilo(Metre)));

  std::int64_t number_of_evaluations = 0;
  FlightPlanOptimizer optimizer(
      flight_plan_,
      FlightPlanOptimizer::LinearCombination(
          {FlightPlanOptimizer::ForCelestialDistance(
               /*celestial=*/&moon,
               /*target_distance=*/2000 * Kilo(Metre)),
           FlightPlanOptimizer::ForΔv()},
          {1, 1e14}),
      [&number_of_evaluations](FlightPlan const&) { ++number_of_evaluations; });

  LOG(INFO) << "Optimizing manœuvre 5";
  auto const manœuvre5 = flight_plan_->GetManœuvre(5);
  EXPECT_OK(optimizer.Optimize(/*index=*/5, 1 * Milli(Metre) / Second));

  EXPECT_THAT(
      manœuvre5.initial_time() - flight_plan_->GetManœuvre(5).initial_time(),
      IsNear(16.8_(1) * Micro(Second)));
  EXPECT_THAT(
      (manœuvre5.Δv() - flight_plan_->GetManœuvre(5).Δv()).Norm(),
      IsNear(1.011_(1) * Metre / Second));

  ComputeFlyby(*flight_plan_, moon, flyby_time, flyby_distance);
  EXPECT_THAT(flyby_time, ResultOf(&TTSecond, "1972-03-27T01:22:33"_DateTime));
  EXPECT_THAT(flyby_distance, IsNear(3339.88_(1) * Kilo(Metre)));
  EXPECT_EQ(146, number_of_evaluations);
}

#if !_DEBUG
struct MetricTestParam {
  MetricTestParam(
      FlightPlanOptimizer::MetricFactory const& metric_factory,
      Matcher<double> const& max_gradient_relative_error,
      Matcher<double> const& max_gateaux_relative_error)
      : metric_factory(metric_factory),
        max_gradient_relative_error(max_gradient_relative_error),
        max_gateaux_relative_error(max_gateaux_relative_error) {}
  FlightPlanOptimizer::MetricFactory const metric_factory;
  Matcher<double> const max_gradient_relative_error;
  Matcher<double> const max_gateaux_relative_error;
};

class MetricTest
    : public ::testing::TestWithParam<MetricTestParam> {
 public:
  using TestNavigationFrame =
      BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation>;

  static Celestial const* earth_celestial() {
    return earth_celestial_.get();
  }

  static std::function<not_null<std::unique_ptr<TestNavigationFrame>>()>
  frame() {
    return []() {
      return make_not_null_unique<TestNavigationFrame>(ephemeris_.get(),
                                                       earth_body_);
    };
  }

 protected:
  MetricTest()
      : epoch_(solar_system_1950_.epoch()),
        navigation_frame_(*frame()()),
        burn_(MakeBurn()) {
    Instant const desired_final_time = epoch_ + 1 * Hour;
    EXPECT_OK(ephemeris_->Prolong(desired_final_time));
    auto const earth_dof = solar_system_1950_.degrees_of_freedom("Earth");
    auto const mars_dof = solar_system_1950_.degrees_of_freedom("Mars");
    auto const midway = Barycentre({earth_dof, mars_dof});
    EXPECT_OK(root_.Append(epoch_, midway));
    flight_plan_ = std::make_unique<FlightPlan>(
        /*initial_mass=*/1 * Kilogram,
        /*initial_time=*/epoch_,
        /*initial_degrees_of_freedom=*/midway,
        desired_final_time,
        ephemeris_.get(),
        Ephemeris<Barycentric>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
            /*max_steps=*/1,
            /*length_integration_tolerance=*/1 * Metre,
            /*speed_integration_tolerance=*/1 * Metre / Second));
    EXPECT_OK(flight_plan_->Insert(burn_, /*index=*/0));
    optimizer_ = std::make_unique<FlightPlanOptimizer>(
        flight_plan_.get(), GetParam().metric_factory);
    metric_ =
        GetParam().metric_factory(optimizer_.get(),
                                  NavigationManœuvre(10 * Kilo(Gram), burn_),
                                  /*index=*/0);
  }

  NavigationManœuvre::Burn MakeBurn() {
    NavigationManœuvre::Intensity intensity;
    intensity.Δv = Velocity<Frenet<Navigation>>(
        {1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second});
    NavigationManœuvre::Timing timing;
    timing.initial_time = epoch_ + 1 * Minute;
    return {intensity,
            timing,
            /*thrust=*/1 * Newton,
            /*specific_impulse=*/1 * Newton * Second / Kilogram,
            make_not_null_unique<TestNavigationFrame>(navigation_frame_),
            /*is_inertially_fixed=*/true};
  }

  static SolarSystem<Barycentric> const solar_system_1950_;
  static not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  static not_null<RotatingBody<Barycentric> const*> const earth_body_;
  static not_null<std::unique_ptr<Celestial>> const earth_celestial_;

  Instant const epoch_;
  TestNavigationFrame const navigation_frame_;
  NavigationManœuvre::Burn const burn_;

  DiscreteTrajectory<Barycentric> root_;
  std::unique_ptr<FlightPlan> flight_plan_;
  std::unique_ptr<FlightPlanOptimizer> optimizer_;
  std::unique_ptr<FlightPlanOptimizer::Metric> metric_;
};

SolarSystem<Barycentric> const MetricTest::solar_system_1950_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2433282_500000000.proto.txt",
    /*ignore_frame=*/true);

not_null<std::unique_ptr<Ephemeris<Barycentric>>> const MetricTest::ephemeris_(
    solar_system_1950_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<Barycentric>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*step=*/10 * Minute)));

not_null<RotatingBody<Barycentric> const*> const MetricTest::earth_body_(
    solar_system_1950_.rotating_body(*ephemeris_, "Earth"));

not_null<std::unique_ptr<Celestial>> const MetricTest::earth_celestial_ = []() {
  auto celestial = std::make_unique<Celestial>(earth_body_);
  celestial->set_trajectory(ephemeris_->trajectory(earth_body_));
  return celestial;
}();

TEST_P(MetricTest, Positive) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate(-100, 100);
  for (int i = 0; i < 100; ++i) {
    FlightPlanOptimizer::HomogeneousArgument const argument(
        std::array{coordinate(random),
                   coordinate(random),
                   coordinate(random),
                   coordinate(random)});
    EXPECT_LE(0, metric_->Evaluate(argument));
  }
}

TEST_P(MetricTest, Gradient) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate(-100, 100);
  std::uniform_real_distribution<double> displacement(-1, 1);
  double max_relative_error = 0.0;
  for (int i = 0; i < 100; ++i) {
    FlightPlanOptimizer::HomogeneousArgument const argument(
        std::array{coordinate(random),
                   coordinate(random),
                   coordinate(random),
                   coordinate(random)});
    for (int j = 0; j < 10; ++j) {
      FlightPlanOptimizer::HomogeneousArgument const Δargument(
          std::array{displacement(random),
                     displacement(random),
                     displacement(random),
                     displacement(random)});
      max_relative_error = std::max(
          max_relative_error,
          RelativeError(
              metric_->Evaluate(argument + Δargument),
              metric_->Evaluate(argument) +
                  TransposedView<FlightPlanOptimizer::HomogeneousArgument>{
                      metric_->EvaluateGradient(argument)} *
                      Δargument));
    }
  }
  EXPECT_THAT(max_relative_error,
              GetParam().max_gradient_relative_error);
}

TEST_P(MetricTest, Gateaux) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate(-100, 100);
  std::uniform_real_distribution<double> displacement(-1, 1);
  double max_relative_error = 0.0;
  for (int i = 0; i < 100; ++i) {
    FlightPlanOptimizer::HomogeneousArgument const argument(
        std::array{coordinate(random),
                   coordinate(random),
                   coordinate(random),
                   coordinate(random)});
    for (int j = 0; j < 10; ++j) {
      FlightPlanOptimizer::HomogeneousArgument const Δargument(
          std::array{displacement(random),
                     displacement(random),
                     displacement(random),
                     displacement(random)});
      max_relative_error = std::max(
          max_relative_error,
          RelativeError(
              metric_->Evaluate(argument + Δargument),
              metric_->Evaluate(argument) +
                  metric_->EvaluateGateauxDerivative(argument, Δargument)));
    }
  }
  EXPECT_THAT(max_relative_error,
              GetParam().max_gateaux_relative_error);
}

INSTANTIATE_TEST_SUITE_P(
    AllMetricTests,
    MetricTest,
    ::testing::Values(
        MetricTestParam(FlightPlanOptimizer::ForCelestialCentre(
                            MetricTest::earth_celestial()),
                        AnyOf(IsNear(1.5e-11_(1)), IsNear(3.9e-11_(1))),
                        AnyOf(IsNear(2.8e-11_(1)), IsNear(7.2e-11_(1)))),
        MetricTestParam(FlightPlanOptimizer::ForCelestialDistance(
                            MetricTest::earth_celestial(),
                            1000 * Kilo(Metre)),
                        AnyOf(IsNear(3.1e-11_(1)), IsNear(7.8e-11_(1))),
                        AnyOf(IsNear(5.6e-11_(1)), IsNear(1.4e-10_(1)))),
        MetricTestParam(FlightPlanOptimizer::ForInclination(
                            MetricTest::earth_celestial(),
                            MetricTest::frame(),
                            45 * Degree),
                        AnyOf(IsNear(1.6e-11_(1)), IsNear(7.8e-11_(1))),
                        AnyOf(IsNear(3.0e-11_(1)), IsNear(1.4e-10_(1)))),
        MetricTestParam(FlightPlanOptimizer::ForΔv(),
                        IsNear(2.6e-3_(1)),
                        IsNear(2.6e-3_(1))),
        MetricTestParam(FlightPlanOptimizer::LinearCombination(
                            {FlightPlanOptimizer::ForCelestialCentre(
                                 MetricTest::earth_celestial()),
                             FlightPlanOptimizer::ForCelestialDistance(
                                 MetricTest::earth_celestial(),
                                 1000 * Kilo(Metre)),
                             FlightPlanOptimizer::ForΔv()},
                            {2, 3, 5}),
                        AnyOf(IsNear(3.1e-11_(1)), IsNear(7.8e-11_(1))),
                        AnyOf(IsNear(5.6e-11_(1)), IsNear(1.4e-10_(1))))));
#endif

}  // namespace ksp_plugin
}  // namespace principia
