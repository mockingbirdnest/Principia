
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/time_scales.hpp"
#include "astronomy/mercury_orbiter.hpp"
#include "base/array.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "base/serialization.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_traject0ry.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/serialization.hpp"
#include "testing_utilities/string_log_sink.hpp"

namespace principia {
namespace interface {

using astronomy::operator""_TT;
using astronomy::MercuryOrbiterInitialDegreesOfFreedom;
using astronomy::MercuryOrbiterInitialTime;
using astronomy::TTSecond;
using astronomy::date_time::DateTime;
using astronomy::date_time::operator""_DateTime;
using base::not_null;
using base::OFStream;
using base::ParseFromBytes;
using base::PullSerializer;
using base::PushDeserializer;
using ksp_plugin::Barycentric;
using ksp_plugin::Plugin;
using physics::DiscreteTraject0ry;
using quantities::Speed;
using quantities::si::Kilo;
using testing_utilities::operator""_⑴;
using testing_utilities::IsNear;
using testing_utilities::ReadFromBinaryFile;
using testing_utilities::ReadLinesFromBase64File;
using testing_utilities::ReadLinesFromHexadecimalFile;
using testing_utilities::StringLogSink;
using testing_utilities::WriteToBinaryFile;
using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::HasSubstr;
using ::testing::Not;
using ::testing::NotNull;
using ::testing::Pair;
using ::testing::SizeIs;
using ::testing::internal::CaptureStderr;
using ::testing::internal::GetCapturedStderr;

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

  // Reads a plugin from a file containing only the "serialized_plugin = "
  // lines, with "serialized_plugin = " dropped.
  static not_null<std::unique_ptr<Plugin const>> ReadPluginFromFile(
      std::filesystem::path const& filename,
      std::string_view const compressor,
      std::string_view const encoder) {
    Plugin const* plugin = nullptr;

    PushDeserializer* deserializer = nullptr;
    auto const lines =
        encoder == "hexadecimal" ? ReadLinesFromHexadecimalFile(filename)
        : encoder == "base64"    ? ReadLinesFromBase64File(filename)
                                  : std::vector<std::string>{};
    CHECK(!lines.empty());

    LOG(ERROR) << "Deserialization starting";
    for (std::string const& line : lines) {
      principia__DeserializePlugin(line.c_str(),
                                   &deserializer,
                                   &plugin,
                                   compressor.data(),
                                   encoder.data());
    }
    principia__DeserializePlugin("",
                                 &deserializer,
                                 &plugin,
                                 compressor.data(),
                                 encoder.data());
    LOG(ERROR) << "Deserialization complete";

    return std::unique_ptr<Plugin const>(plugin);
  }

  // Writes a plugin to a file.
  static void WritePluginToFile(
      std::filesystem::path const& filename,
      std::string_view const compressor,
      std::string_view const encoder,
      not_null<std::unique_ptr<Plugin const>> plugin) {
    OFStream file(filename);
    PullSerializer* serializer = nullptr;
    char const* b64 = nullptr;

    LOG(ERROR) << "Serialization starting";
    for (;;) {
      b64 = principia__SerializePlugin(plugin.get(),
                                       &serializer,
                                       preferred_compressor,
                                       preferred_encoder);
      if (b64 == nullptr) {
        break;
      }
      file << b64 << "\n";
      principia__DeleteString(&b64);
    }
    LOG(ERROR) << "Serialization complete";

    Plugin const* released_plugin = plugin.release();
    principia__DeletePlugin(&released_plugin);
  }

  static void WriteAndReadBack(
      not_null<std::unique_ptr<Plugin const>> plugin1) {
    // Write the plugin to a new file with the preferred format.
    WritePluginToFile(TEMP_DIR / "serialized_plugin.proto.b64",
                      preferred_compressor,
                      preferred_encoder,
                      std::move(plugin1));

    // Read the plugin from the new file to make sure that it's fine.
    auto plugin2 = ReadPluginFromFile(TEMP_DIR / "serialized_plugin.proto.b64",
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
  EXPECT_THAT(TTSecond(test->history()->front().time),
              Eq("1970-08-14T08:03:18"_DateTime));
  EXPECT_THAT(TTSecond(test->psychohistory()->back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  EXPECT_FALSE(test->has_flight_plan());

  auto const ifnity = plugin->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");
  EXPECT_THAT(ifnity->name(), Eq("IFNITY-5.2"));
  EXPECT_THAT(TTSecond(ifnity->history()->front().time),
              Eq("1970-08-14T08:03:46"_DateTime));
  EXPECT_THAT(TTSecond(ifnity->psychohistory()->back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  ASSERT_TRUE(ifnity->has_flight_plan());
  EXPECT_THAT(ifnity->flight_plan().number_of_manœuvres(), Eq(16));
  std::vector<std::pair<DateTime, Speed>> manœuvre_ignition_tt_seconds_and_Δvs;
  for (int i = 0; i < ifnity->flight_plan().number_of_manœuvres(); ++i) {
    manœuvre_ignition_tt_seconds_and_Δvs.emplace_back(
        TTSecond(ifnity->flight_plan().GetManœuvre(i).initial_time()),
        ifnity->flight_plan().GetManœuvre(i).Δv().Norm());
  }
  // The flight plan only covers the inner solar system (this is probably
  // because of #3035).
  // It also differs from https://youtu.be/7BDxZV7UD9I?t=439.
  // TODO(egg): Compute the flybys and figure out what exactly is going on in
  // this flight plan.
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
              AllOf(HasSubstr("pre-Haar"), Not(HasSubstr(u8"pre-Gröbner"))));
  auto const& orbiter =
      *plugin->GetVessel("e180ca12-492f-45bf-a194-4c5255aec8a0");
  EXPECT_THAT(orbiter.name(), Eq("Mercury Orbiter 1"));
  auto const begin = orbiter.history()->begin();
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
              IsNear(176'400'999_⑴ * Kilo(Metre)));
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
              IsNear(19'163_⑴ * Kilo(Metre)));

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
              AllOf(HasSubstr(u8"pre-Ζήνων"), Not(HasSubstr("pre-Haar"))));

  // The vessel with the longest history.
  auto const& vessel =
      *plugin->GetVessel("77ddea45-47ee-48c0-aee9-d55cdb35ffcd");
  auto history = vessel.history();
  auto psychohistory = vessel.psychohistory();
  EXPECT_THAT(*history, SizeIs(435'927));
  EXPECT_THAT(*psychohistory, SizeIs(3));

  // Evaluate a point in each of the two segments.
  EXPECT_THAT(history->EvaluateDegreesOfFreedom("1957-10-04T19:28:34"_TT),
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
    auto const serialized_message = base::SerializeAsBytes(message);
    WriteToBinaryFile(TEMP_DIR / "trajectory_3136.proto.bin",
                      serialized_message.get());
  }

  // Deserialize the temporary file to make sure that it's valid.
  {
    auto const serialized_message =
        ReadFromBinaryFile(TEMP_DIR / "trajectory_3136.proto.bin");
    auto const message =
        ParseFromBytes<serialization::DiscreteTrajectory>(serialized_message);
    auto const trajectory = DiscreteTraject0ry<Barycentric>::ReadFromMessage(
        message, /*tracked=*/{&history, &psychohistory});
    EXPECT_THAT(*history, SizeIs(435'927));
    EXPECT_THAT(*psychohistory, SizeIs(3));
  }

  // Make sure that we can upgrade, save, and reload.
  WriteAndReadBack(std::move(plugin));
}

// Use for debugging saves given by users.
TEST_F(PluginCompatibilityTest, DISABLED_SECULAR_Debug) {
  CheckSaveCompatibility(
      R"(P:\Public Mockingbird\Principia\Saves\2685\five-minute-scene-change-neptune.txt)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
}

}  // namespace interface
}  // namespace principia
