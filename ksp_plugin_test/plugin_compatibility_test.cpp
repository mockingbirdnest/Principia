#include "ksp_plugin/plugin.hpp"

#include <map>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/date_time.hpp"
#include "astronomy/mercury_orbiter.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "base/serialization.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/interface.hpp"  // 🧙 For interface functions.
#include "ksp_plugin_test/plugin_io.hpp"
#include "numerics/fma.hpp"
#include "physics/apsides.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/serialization.hpp"
#include "testing_utilities/string_log_sink.hpp"

namespace principia {
namespace interface {

using ::testing::AllOf;
using ::testing::Bool;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::HasSubstr;
using ::testing::IsEmpty;
using ::testing::Not;
using ::testing::NotNull;
using ::testing::Pair;
using ::testing::ResultOf;
using ::testing::SizeIs;
using ::testing::internal::CaptureStderr;
using ::testing::internal::GetCapturedStderr;
using namespace principia::astronomy::_date_time;
using namespace principia::astronomy::_mercury_orbiter;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::base::_serialization;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin_test::_plugin_io;
using namespace principia::numerics::_fma;
using namespace principia::physics::_apsides;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_serialization;
using namespace principia::testing_utilities::_string_log_sink;

using namespace std::chrono_literals;

const char preferred_compressor[] = "gipfeli";
const char preferred_encoder[] = "base64";

class PluginCompatibilityTest : public testing::Test {
 protected:
  PluginCompatibilityTest()
      : stderrthreshold_(FLAGS_stderrthreshold) {
    google::SetStderrLogging(google::WARNING);
  }

  ~PluginCompatibilityTest() override {
    google::SetStderrLogging(stderrthreshold_);
  }

  static void WriteAndReadBack(
      not_null<std::unique_ptr<Plugin const>> plugin1) {
    // There may be solidi in the path due to parameterized tests, so we remove
    // them.
    std::string const sanitized_name = absl::StrCat(
        absl::StrReplaceAll(
            testing::UnitTest::GetInstance()->current_test_info()->name(),
            {{"/", "__"}}),
        "_serialized_plugin.proto.b64");
    // Write the plugin to a new file with the preferred format.
    WritePluginToFile(TEMP_DIR / sanitized_name,
                      preferred_compressor,
                      preferred_encoder,
                      std::move(plugin1));

    // Read the plugin from the new file to make sure that it's fine.
    auto plugin2 = ReadPluginFromFile(TEMP_DIR / sanitized_name,
                                      preferred_compressor,
                                      preferred_encoder);
  }

  static void CheckSaveCompatibility(std::filesystem::path const& filename,
                                     std::string_view const compressor,
                                     std::string_view const encoder) {
    // Read a plugin from the given file.
    auto plugin = ReadPluginFromFile(filename, compressor, encoder);

    WriteAndReadBack(std::move(plugin));
  }

  int const stderrthreshold_;
};

TEST_F(PluginCompatibilityTest, PreCartan) {
  // This space for rent.
}

TEST_F(PluginCompatibilityTest, PreCohen) {
  StringLogSink log_warning(google::WARNING);
  CheckSaveCompatibility(
      SOLUTION_DIR / "ksp_plugin_test" / "saves" / "3039.proto.hex",
      /*compressor=*/"",
      /*decoder=*/"hexadecimal");
  EXPECT_THAT(
      log_warning.string(),
      AllOf(
          HasSubstr(
              "pre-Cohen ContinuousTrajectory"),  // Regression test for #3039.
          HasSubstr("pre-Cauchy"),                // The save is even older.
          Not(HasSubstr("pre-Cartan"))));         // But not *that* old.
}

#if !_DEBUG

TEST_F(PluginCompatibilityTest, Reach) {
  StringLogSink log_warning(google::WARNING);
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      SOLUTION_DIR / "ksp_plugin_test" / "saves" / "3072.proto.b64",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
  EXPECT_THAT(log_warning.string(),
              AllOf(HasSubstr("pre-Galileo"), Not(HasSubstr("pre-Frobenius"))));
  auto const test = plugin->GetVessel("f2d77873-4776-4809-9dfb-de9e7a0620a6");
  EXPECT_THAT(test->name(), Eq("TEST"));
  EXPECT_THAT(TTSecond(test->trajectory().front().time),
              Eq("1970-08-14T08:03:18"_DateTime));
  EXPECT_THAT(TTSecond(test->psychohistory()->back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  EXPECT_FALSE(test->has_flight_plan());

  auto const ifnity = plugin->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");
  EXPECT_THAT(ifnity->name(), Eq("IFNITY-5.2"));
  EXPECT_THAT(TTSecond(ifnity->trajectory().front().time),
              Eq("1970-08-14T08:03:46"_DateTime));
  EXPECT_THAT(TTSecond(ifnity->psychohistory()->back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  ASSERT_TRUE(ifnity->has_flight_plan());
  ifnity->ReadFlightPlanFromMessage();
  FlightPlan const& flight_plan = ifnity->flight_plan();
  EXPECT_THAT(flight_plan.desired_final_time(),
              ResultOf(&TTSecond, Eq("1980-03-14T09:37:01"_DateTime)));
  EXPECT_THAT(flight_plan.actual_final_time(),
              ResultOf(&TTSecond, Eq("1980-03-14T09:37:01"_DateTime)));
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
  // The flight plan only covers the inner solar system (this is probably
  // because of #3035).
  // The manœuvres differ from those in https://youtu.be/7BDxZV7UD9I?t=439.
  EXPECT_THAT(manœuvre_ignition_tt_seconds_and_Δvs,
              ElementsAre(Pair("1970-08-14T09:34:49"_DateTime,
                               3.80488671073918022e+03 * (Metre / Second)),
                          Pair("1970-08-15T13:59:24"_DateTime,
                               3.04867185471741759e-04 * (Metre / Second)),
                          Pair("1970-12-22T07:48:21"_DateTime,
                               1.58521291818444873e-03 * (Metre / Second)),
                          Pair("1971-01-08T17:36:55"_DateTime,
                               1.40000000034068623e-03 * (Metre / Second)),
                          Pair("1971-07-02T17:16:00"_DateTime,
                               1.00000000431022681e-04 * (Metre / Second)),
                          Pair("1971-09-06T03:27:33"_DateTime,
                               1.78421858738381537e-03 * (Metre / Second)),
                          Pair("1972-02-13T22:47:26"_DateTime,
                               7.72606625794511597e-04 * (Metre / Second)),
                          Pair("1972-03-25T16:30:19"_DateTime,
                               5.32846131747503372e-03 * (Metre / Second)),
                          Pair("1972-12-24T04:09:32"_DateTime,
                               3.45000000046532824e-03 * (Metre / Second)),
                          Pair("1973-06-04T01:59:07"_DateTime,
                               9.10695453328359134e-03 * (Metre / Second)),
                          Pair("1973-07-09T06:07:17"_DateTime,
                               4.49510921430966881e-01 * (Metre / Second)),
                          Pair("1973-09-10T03:59:44"_DateTime,
                               1.00000000431022681e-04 * (Metre / Second)),
                          Pair("1974-11-20T17:34:27"_DateTime,
                               5.10549409572428781e-01 * (Metre / Second)),
                          Pair("1975-10-07T01:29:45"_DateTime,
                               2.86686518692948443e-02 * (Metre / Second)),
                          Pair("1975-12-29T21:27:13"_DateTime,
                               1.00404183285598275e-03 * (Metre / Second)),
                          Pair("1977-07-28T22:47:53"_DateTime,
                               1.39666705839172456e-01 * (Metre / Second))));

  // Compute the flybys.
  std::map<Instant, std::string> flyby_map;
  for (Index index = 0; plugin->HasCelestial(index); ++index) {
    Celestial const& celestial = plugin->GetCelestial(index);
    auto const& celestial_trajectory = celestial.trajectory();
    auto const& flight_plan_trajectory = flight_plan.GetAllSegments();
    DiscreteTrajectory<Barycentric> apoapsides;
    DiscreteTrajectory<Barycentric> periapsides;

    // The begin time avoid spurious periapsides right after the launch.
    ComputeApsides(celestial_trajectory,
                   flight_plan_trajectory,
                   flight_plan_trajectory.upper_bound("1970-08-15T00:00:00"_TT),
                   flight_plan_trajectory.end(),
                   /*t_max=*/InfiniteFuture,
                   /*max_points=*/100,
                   apoapsides,
                   periapsides);
    auto const radius = celestial.body()->mean_radius();
    for (auto const& [time, _] : periapsides) {
      if ((celestial_trajectory.EvaluatePosition(time) -
           flight_plan_trajectory.EvaluatePosition(time))
              .Norm() < 50 * radius) {
        flyby_map[time] = celestial.body()->name();
      }
    }
  }
  std::vector<std::pair<Instant, std::string>> flybys(flyby_map.begin(),
                                                      flyby_map.end());

#if PRINCIPIA_COMPILER_MSVC
  EXPECT_THAT(
      flybys,
      ElementsAre(
          Pair(ResultOf(&TTSecond, "1970-12-23T07:16:42"_DateTime), "Venus"),
          Pair(ResultOf(&TTSecond, "1971-08-29T23:33:54"_DateTime), "Mars"),
          Pair(ResultOf(&TTSecond, "1972-03-26T11:23:30"_DateTime), "Earth"),
          Pair(ResultOf(&TTSecond, "1972-03-27T01:02:40"_DateTime), "Moon"),
          Pair(ResultOf(&TTSecond, "1972-11-02T21:15:50"_DateTime), "Venus"),
          Pair(ResultOf(&TTSecond, "1973-06-15T13:17:25"_DateTime), "Venus"),
          Pair(ResultOf(&TTSecond, "1974-07-25T22:45:52"_DateTime), "Mercury"),
          Pair(ResultOf(&TTSecond, "1974-09-05T08:27:45"_DateTime), "Venus"),
          Pair(ResultOf(&TTSecond, "1975-04-18T00:42:26"_DateTime), "Venus"),
          Pair(ResultOf(&TTSecond, "1976-04-26T17:36:24"_DateTime), "Venus"),
          Pair(  // The video has             21:57.
              ResultOf(&TTSecond, "1978-08-07T21:58:49"_DateTime),
              "Earth"),
          Pair(  // The video has             07:52.
              ResultOf(&TTSecond, "1980-02-17T22:18:43"_DateTime),
              "Jupiter")));
#elif OS_MACOSX
  EXPECT_THAT(
      flybys,
      ElementsAre(
          Pair(ResultOf(&TTSecond, "1970-12-23T07:16:42"_DateTime), "Venus"),
          Pair(ResultOf(&TTSecond, "1971-08-29T23:34:21"_DateTime), "Mars")));
#elif OS_LINUX
  EXPECT_THAT(
      flybys,
      ElementsAre(
          Pair(ResultOf(&TTSecond, "1970-12-23T07:16:42"_DateTime), "Venus"),
          Pair(ResultOf(&TTSecond, "1971-08-29T23:34:20"_DateTime), "Mars")));
#else
#error Unknown OS
#endif

  // Make sure that we can upgrade, save, and reload.
  WriteAndReadBack(std::move(plugin));
}
#endif

TEST_F(PluginCompatibilityTest, DISABLED_Butcher) {
  StringLogSink log_warning(google::WARNING);
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      R"(P:\Public Mockingbird\Principia\Saves\1119\1119.proto.b64)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
  EXPECT_THAT(log_warning.string(),
              AllOf(HasSubstr("pre-Haar"), Not(HasSubstr("pre-Gröbner"))));
  auto const& orbiter =
      *plugin->GetVessel("e180ca12-492f-45bf-a194-4c5255aec8a0");
  EXPECT_THAT(orbiter.name(), Eq("Mercury Orbiter 1"));
  auto const begin = orbiter.trajectory().begin();
  EXPECT_THAT(begin->time,
              Eq("1966-05-10T00:14:03"_TT + 0.0879862308502197 * Second));
  EXPECT_THAT(begin->degrees_of_freedom,
              Eq(DegreesOfFreedom<Barycentric>(
                  Barycentric::origin + Displacement<Barycentric>(
                                            {-9.83735958466250000e+10 * Metre,
                                             -1.05659916408781250e+11 * Metre,
                                             -4.58171358797500000e+10 * Metre}),
                  Velocity<Barycentric>(
                      {+2.18567382812500000e+04 * (Metre / Second),
                       -1.76616533203125000e+04 * (Metre / Second),
                       -7.76112133789062500e+03 * (Metre / Second)}))));

  auto const& mercury = plugin->GetCelestial(2);
  EXPECT_THAT(mercury.body()->name(), Eq("Mercury"));

  plugin->RequestReanimation(begin->time);

  while (mercury.trajectory().t_min() > begin->time) {
    absl::SleepFor(absl::Milliseconds(1));
  }

  // The history goes back far enough that we are still on our way to Mercury at
  // the beginning.
  EXPECT_THAT((begin->degrees_of_freedom.position() -
               mercury.trajectory().EvaluatePosition(begin->time)).Norm(),
              IsNear(176'400'999_(1) * Kilo(Metre)));
  EXPECT_THAT(begin->time,
              Eq("1966-05-10T00:14:03"_TT + 0.0879862308502197 * Second));
  EXPECT_THAT(begin->degrees_of_freedom,
              Eq(DegreesOfFreedom<Barycentric>(
                  Barycentric::origin + Displacement<Barycentric>(
                                            {-9.83735958466250000e+10 * Metre,
                                             -1.05659916408781250e+11 * Metre,
                                             -4.58171358797500000e+10 * Metre}),
                  Velocity<Barycentric>(
                      {+2.18567382812500000e+04 * (Metre / Second),
                       -1.76616533203125000e+04 * (Metre / Second),
                       -7.76112133789062500e+03 * (Metre / Second)}))));

  // We arrive in late August.  Check the state in the beginning of September.
  auto const it = orbiter.trajectory().lower_bound("1966-09-01T00:00:00"_TT);
  EXPECT_THAT(it->time, Eq(MercuryOrbiterInitialTime));
  EXPECT_THAT(it->degrees_of_freedom,
              Eq(MercuryOrbiterInitialDegreesOfFreedom<Barycentric>));
  EXPECT_THAT((it->degrees_of_freedom.position() -
               mercury.trajectory().EvaluatePosition(it->time)).Norm(),
              IsNear(19'163_(1) * Kilo(Metre)));

  // Make sure that we can upgrade, save, and reload.
  WriteAndReadBack(std::move(plugin));
}

TEST_F(PluginCompatibilityTest, DISABLED_Lpg) {
  StringLogSink log_warning(google::WARNING);
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      R"(P:\Public Mockingbird\Principia\Saves\3136\3136.proto.b64)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
  EXPECT_THAT(log_warning.string(),
              AllOf(HasSubstr("pre-Hamilton"), Not(HasSubstr("pre-Haar"))));

  // The vessel with the longest history.
  auto const& vessel =
      *plugin->GetVessel("77ddea45-47ee-48c0-aee9-d55cdb35ffcd");
  auto const& trajectory = vessel.trajectory();
  auto history = trajectory.segments().begin();
  auto psychohistory = vessel.psychohistory();
  EXPECT_THAT(trajectory, SizeIs(435'929));
  EXPECT_THAT(*history, SizeIs(435'927));
  EXPECT_THAT(*psychohistory, SizeIs(3));

  // Evaluate a point in each of the two segments.
  EXPECT_THAT(trajectory.EvaluateDegreesOfFreedom("1957-10-04T19:28:34"_TT),
              Eq(DegreesOfFreedom<Barycentric>(
                  Barycentric::origin + Displacement<Barycentric>(
                                            {+1.47513683827317657e+11 * Metre,
                                             +2.88696086355042419e+10 * Metre,
                                             +1.24740082262952404e+10 * Metre}),
                  Velocity<Barycentric>(
                      {-6.28845231836519179e+03 * (Metre / Second),
                       +2.34046542233168329e+04 * (Metre / Second),
                       +4.64410011408655919e+03 * (Metre / Second)}))));
  EXPECT_THAT(psychohistory->EvaluateDegreesOfFreedom("1958-10-07T09:38:30"_TT),
              Eq(DegreesOfFreedom<Barycentric>(
                  Barycentric::origin + Displacement<Barycentric>(
                                            {+1.45814173315801941e+11 * Metre,
                                             +3.45409490426372147e+10 * Metre,
                                             +1.49445864962450924e+10 * Metre}),
                  Velocity<Barycentric>(
                      {-8.70708379504568074e+03 * (Metre / Second),
                       +2.61488327506437054e+04 * (Metre / Second),
                       +1.90319283138508908e+04 * (Metre / Second)}))));

  // Serialize the history and psychohistory to a temporary file.
  {
    serialization::DiscreteTrajectory message;
    vessel.trajectory().WriteToMessage(
        &message, /*tracked=*/{history, psychohistory}, /*exact=*/{});
    auto const serialized_message = SerializeAsBytes(message);
    WriteToBinaryFile(TEMP_DIR / "trajectory_3136.proto.bin",
                      serialized_message.get());
  }

  // Deserialize the temporary file to make sure that it's valid.
  {
    auto const serialized_message =
        ReadFromBinaryFile(TEMP_DIR / "trajectory_3136.proto.bin");
    auto const message =
        ParseFromBytes<serialization::DiscreteTrajectory>(serialized_message);
    auto const trajectory = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message, /*tracked=*/{&history, &psychohistory});
    EXPECT_THAT(trajectory, SizeIs(435'929));
    EXPECT_THAT(*history, SizeIs(435'927));
    EXPECT_THAT(*psychohistory, SizeIs(3));
  }

  // Make sure that we can upgrade, save, and reload.
  WriteAndReadBack(std::move(plugin));
}

TEST_F(PluginCompatibilityTest, DISABLED_Egg) {
  StringLogSink log_warning(google::WARNING);
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      R"(P:\Public Mockingbird\Principia\Saves\3136\3136b.proto.b64)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
  EXPECT_THAT(log_warning.string(),
              AllOf(HasSubstr("pre-Hamilton"), Not(HasSubstr("pre-Haar"))));

  auto& mutable_plugin = const_cast<Plugin&>(*plugin);

  // This would fail if segment iterators were invalidated during part
  // deserialization.
  mutable_plugin.AdvanceTime(
      mutable_plugin.GameEpoch() + 133218.91123694609 * Second,
      295.52698460805016 * Degree);
  mutable_plugin.CatchUpVessel("1e07aaa2-d1f8-4f6d-8b32-495b46109d98");

  // Make sure that we can upgrade, save, and reload.
  WriteAndReadBack(std::move(plugin));
}

TEST_F(PluginCompatibilityTest, PreHardy) {
  StringLogSink log_warning(google::WARNING);
  CheckSaveCompatibility(
      SOLUTION_DIR / "ksp_plugin_test" / "saves" / "3244.proto.b64",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
  // Regression test for #3244.
  EXPECT_THAT(log_warning.string(),
              AllOf(HasSubstr("pre-Hardy DiscreteTrajectorySegment"),
                    Not(HasSubstr("pre-Hamilton"))));
}

#if !_DEBUG
TEST_F(PluginCompatibilityTest, 3273) {
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      SOLUTION_DIR / "ksp_plugin_test" / "saves" / "3273.proto.b64",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
  auto const& vessel =
      plugin->GetVessel("ae3aa35c-f33c-486e-a59c-aee43954dc30");
  vessel->ReadFlightPlanFromMessage();
  EXPECT_THAT(vessel->flight_plan().number_of_manœuvres(), Eq(2));
  while (vessel->flight_plan().analysis(2) == nullptr) {
    LOG(ERROR) << static_cast<int>(
                      vessel->flight_plan().progress_of_analysis(2) * 100)
               << "%";
    std::this_thread::sleep_for(1s);
  }
  auto const& analysis = *vessel->flight_plan().analysis(2);
  auto const Ωʹ = analysis.elements()->nodal_precession();
  auto const Ωʹᴛ = analysis.primary()->angular_frequency();
  auto const nd = 2 * π * Radian / analysis.elements()->nodal_period();
  double const κ = nd / (Ωʹᴛ - Ωʹ);
  EXPECT_THAT(vessel->flight_plan().analysis(2)->recurrence(),
              Eq(std::nullopt));
}
#endif

// Use for debugging saves given by users.
TEST_F(PluginCompatibilityTest, DISABLED_SECULAR_Debug) {
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      R"(P:\Public Mockingbird\Principia\Saves\3203\wip.proto.b64)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
}

}  // namespace interface
}  // namespace principia
