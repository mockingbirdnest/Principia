
#include <filesystem>
#include <fstream>
#include <map>
#include <string>

#include "astronomy/time_scales.hpp"
#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "google/protobuf/io/coded_stream.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/apsides.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/serialization.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using astronomy::TTSecond;
using astronomy::operator""_TT;
using astronomy::date_time::DateTime;
using astronomy::date_time::operator""_DateTime;
using base::Array;
using base::HexadecimalEncoder;
using base::UniqueArray;
using geometry::Bivector;
using geometry::Trivector;
using geometry::Vector;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::ComputeApsides;
using quantities::Length;
using quantities::Speed;
using quantities::si::Hour;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using ::testing::AllOf;
using ::testing::AnyOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::IsFalse;
using ::testing::IsTrue;
using ::testing::Lt;

class PluginCompatibilityTest : public testing::Test {
protected:
  PluginCompatibilityTest() {
    google::LogToStderr();
  }

  std::unique_ptr<Plugin const> ReadPluginFromBase64File(
      std::filesystem::path const& filename) {
    LOG(INFO) << "Reading file " << filename << u8"…";
    auto const lines = testing_utilities::ReadLinesFromBase64File(filename);
    base::PushDeserializer* deserializer = nullptr;
    Plugin const* plugin;
    LOG(INFO) << lines.size() << u8" lines. Deserializing…";
    for (std::string const& line : lines) {
      interface::principia__DeserializePlugin(line.c_str(),
                                              &deserializer,
                                              &plugin,
                                              /*compressor=*/"gipfeli",
                                              "base64");
    }
    interface::principia__DeserializePlugin("",
                                            &deserializer,
                                            &plugin,
                                            /*compressor=*/"gipfeli",
                                            "base64");
    LOG(INFO) << "Deserialization complete.";
    return std::unique_ptr<Plugin const>(plugin);
  }
};

TEST_F(PluginCompatibilityTest, PreCartan) {
  // This space for rent.
}

TEST_F(PluginCompatibilityTest, Reach) {
  std::unique_ptr<Plugin const> const plugin =
      ReadPluginFromBase64File(SOLUTION_DIR / "ksp_plugin_test" / "Reach.sfs");
  auto const test = plugin->GetVessel("f2d77873-4776-4809-9dfb-de9e7a0620a6");
  EXPECT_THAT(test->name(), Eq("TEST"));
  EXPECT_THAT(TTSecond(test->psychohistory().front().time),
              Eq("1970-08-14T08:03:18"_DateTime));
  EXPECT_THAT(TTSecond(test->psychohistory().back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  EXPECT_THAT(test->has_flight_plan(), IsFalse());

  auto const infnity =
      plugin->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");
  EXPECT_THAT(infnity->name(), Eq("IFNITY-5.2"));
  EXPECT_THAT(TTSecond(infnity->psychohistory().front().time),
              Eq("1970-08-14T08:03:46"_DateTime));
  EXPECT_THAT(TTSecond(infnity->psychohistory().back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  EXPECT_THAT(infnity->has_flight_plan(), IsTrue());
  EXPECT_THAT(infnity->flight_plan().number_of_manœuvres(), Eq(16));
  std::vector<DateTime> manœuvre_tt_seconds;
  std::vector<Speed> manœuvre_Δvs;
  for (int i = 0; i < infnity->flight_plan().number_of_manœuvres(); ++i) {
    manœuvre_tt_seconds.push_back(
        TTSecond(infnity->flight_plan().GetManœuvre(i).initial_time()));
    manœuvre_Δvs.push_back(infnity->flight_plan().GetManœuvre(i).Δv().Norm());
  }
  EXPECT_THAT(manœuvre_tt_seconds,
              ElementsAre("1970-08-14T09:34:49"_DateTime,
                          "1970-08-15T13:59:24"_DateTime,
                          "1970-12-22T07:48:21"_DateTime,
                          "1971-01-08T17:36:55"_DateTime,
                          "1971-07-02T17:16:00"_DateTime,
                          "1971-09-06T03:27:33"_DateTime,
                          "1972-02-13T22:47:26"_DateTime,
                          "1972-03-25T16:30:19"_DateTime,
                          "1972-12-24T04:09:32"_DateTime,
                          "1973-06-04T01:59:07"_DateTime,
                          "1973-07-09T06:07:17"_DateTime,
                          "1973-09-10T03:59:44"_DateTime,
                          "1974-11-20T17:34:27"_DateTime,
                          "1975-10-07T01:29:45"_DateTime,
                          "1975-12-29T21:27:13"_DateTime,
                          "1977-07-28T22:47:53"_DateTime));
  EXPECT_THAT(manœuvre_Δvs,
              ElementsAre(+3.80488671073918022e+03 * (Metre / Second),
                          +3.04867185471741759e-04 * (Metre / Second),
                          +1.58521291818444873e-03 * (Metre / Second),
                          +1.40000000034068623e-03 * (Metre / Second),
                          +1.00000000431022681e-04 * (Metre / Second),
                          +1.78421858738381537e-03 * (Metre / Second),
                          +7.72606625794511597e-04 * (Metre / Second),
                          +5.32846131747503372e-03 * (Metre / Second),
                          +3.45000000046532824e-03 * (Metre / Second),
                          +9.10695453328359134e-03 * (Metre / Second),
                          +4.49510921430966881e-01 * (Metre / Second),
                          +1.00000000431022681e-04 * (Metre / Second),
                          +5.10549409572428781e-01 * (Metre / Second),
                          +2.86686518692948443e-02 * (Metre / Second),
                          +1.00404183285598275e-03 * (Metre / Second),
                          +1.39666705839172456e-01 * (Metre / Second)));
  for (Vessel const* vessel : {test, infnity}) {
    LOG(ERROR) << vessel->name() << ":";
    if (vessel->has_flight_plan()) {
      auto& flight_plan = vessel->flight_plan();
      LOG(ERROR) << flight_plan.number_of_manœuvres() << u8" manœuvres";
      LOG(ERROR) << "Flight plan: " << TTSecond(flight_plan.initial_time())
                 << " .. " << TTSecond(flight_plan.actual_final_time());
      auto adaptive_step_parameters = flight_plan.adaptive_step_parameters();
      adaptive_step_parameters.set_max_steps(
          std::numeric_limits<int64_t>::max());
      flight_plan.SetAdaptiveStepParameters(
          adaptive_step_parameters,
          flight_plan.generalized_adaptive_step_parameters());
      for (;;) {
        if (flight_plan.SetDesiredFinalTime("1989-07-14T12:00:00"_TT).ok()) {
          break;
        }
        LOG(ERROR) << flight_plan.actual_final_time();
        LOG(ERROR) << "Extended to "
                   << TTSecond(flight_plan.actual_final_time());
      }
      LOG(ERROR) << "Extended to " << TTSecond(flight_plan.actual_final_time());
      for (int i = 0; i < flight_plan.number_of_manœuvres(); ++i) {
        Instant const t = flight_plan.GetManœuvre(i).initial_time();
        LOG(ERROR) << flight_plan.GetManœuvre(i).Δv().Norm() << " at "
                   << TTSecond(t)
                   << t - astronomy::internal_time_scales::DateTimeAsTT(
                              TTSecond(t));
        serialization::DynamicFrame frame;
        flight_plan.GetManœuvre(i).frame()->WriteToMessage(&frame);
        auto const& manœuvre_frame = dynamic_cast<
            BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation> const&>(
            *flight_plan.GetManœuvre(i).frame());
        LOG(ERROR) << manœuvre_frame.centre()->name();
        DiscreteTrajectory<Barycentric>::Iterator begin;
        DiscreteTrajectory<Barycentric>::Iterator end;
        // Loop over the preceding coast, the current burn, and the final coast
        // if this is the last manœuvre.
        for (int const j : (i == flight_plan.number_of_manœuvres() - 1
                                ? std::vector{0, 1, 2}
                                : std::vector{0, 1})) {
          flight_plan.GetSegment(2 * i + j, begin, end);
          for (int i = 0; i < 100; ++i) {
            if (plugin->HasCelestial(i)) {
              auto const& celestial = plugin->GetCelestial(i);
              DiscreteTrajectory<Barycentric> apoapsides;
              DiscreteTrajectory<Barycentric> periapsides;
              ComputeApsides(celestial.trajectory(),
                             begin,
                             end,
                             /*max_points=*/std::numeric_limits<int>::max(),
                             apoapsides,
                             periapsides);
              for (auto const periapsis : periapsides) {
                auto const distance =
                    (celestial.trajectory().EvaluatePosition(periapsis.time) -
                     periapsis.degrees_of_freedom.position())
                        .Norm();
                std::set<std::string_view> const gas_giants{
                    "Jupiter", "Saturn", "Uranus", "Neptune"};
                Length const threshold =
                    (gas_giants.contains(celestial.body()->name()) ||
                     (celestial.parent() != nullptr &&
                      gas_giants.contains(celestial.parent()->body()->name())))
                        ? 1e7 * Kilo(Metre)
                        : 1e5 * Kilo(Metre);
                if (distance < threshold) {
                  LOG(ERROR) << TTSecond(periapsis.time) << ": "
                             << distance / Kilo(Metre) << " km from "
                             << celestial.body()->name();
                }
              }
            }
          }
        }
      }
    } else {
      LOG(ERROR) << "No flight plan.";
    }
    LOG(ERROR) << "Psychohistory range: "
               << TTSecond(vessel->psychohistory().front().time) << " .. "
               << TTSecond(vessel->psychohistory().back().time);
  }
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
