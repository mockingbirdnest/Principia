
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/time_scales.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "testing_utilities/serialization.hpp"

namespace principia {
namespace interface {

using astronomy::operator""_TT;
using astronomy::TTSecond;
using astronomy::date_time::DateTime;
using astronomy::date_time::operator""_DateTime;
using base::not_null;
using base::OFStream;
using base::PullSerializer;
using base::PushDeserializer;
using ksp_plugin::Plugin;
using quantities::Speed;
using testing_utilities::ReadLinesFromBase64File;
using testing_utilities::ReadLinesFromHexadecimalFile;
using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::HasSubstr;
using ::testing::Not;
using ::testing::NotNull;
using ::testing::Pair;
using ::testing::internal::CaptureStderr;
using ::testing::internal::GetCapturedStderr;

const char preferred_compressor[] = "gipfeli";
const char preferred_encoder[] = "base64";

class StringLogSink : google::LogSink {
 public:
  explicit StringLogSink(google::LogSeverity const minimal_severity)
      : minimal_severity_(minimal_severity) {
    google::AddLogSink(this);
  }

  ~StringLogSink() {
    google::RemoveLogSink(this);
  }

  void send(google::LogSeverity const severity,
            char const* const full_filename,
            char const* const base_filename,
            int const line,
            tm const* const tm_time,
            const char* const message,
            size_t const message_len) override {
    if (severity < minimal_severity_) {
      return;
    }
    absl::MutexLock lock(&mutex_);
    absl::StrAppend(
        &string_,
        ToString(severity, base_filename, line, tm_time, message, message_len));
  }

  std::string& string() {
    return string_;
  }

 private:
  google::LogSeverity const minimal_severity_;
  absl::Mutex mutex_;
  std::string string_ GUARDED_BY(mutex_);
};

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

  int stderrthreshold_;
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
  EXPECT_THAT(TTSecond(test->psychohistory().front().time),
              Eq("1970-08-14T08:03:18"_DateTime));
  EXPECT_THAT(TTSecond(test->psychohistory().back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  EXPECT_FALSE(test->has_flight_plan());

  auto const ifnity = plugin->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");
  EXPECT_THAT(ifnity->name(), Eq("IFNITY-5.2"));
  EXPECT_THAT(TTSecond(ifnity->psychohistory().front().time),
              Eq("1970-08-14T08:03:46"_DateTime));
  EXPECT_THAT(TTSecond(ifnity->psychohistory().back().time),
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

// Use for debugging saves given by users.
TEST_F(PluginCompatibilityTest, DISABLED_SECULAR_Debug) {
  CheckSaveCompatibility(
      R"(P:\Public Mockingbird\Principia\Saves\2685\five-minute-scene-change-neptune.txt)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
}

}  // namespace interface
}  // namespace principia
