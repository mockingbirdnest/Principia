
#include "journal/player.hpp"

#include <list>
#include <string>
#include <thread>
#include <vector>

#include "benchmark/benchmark.h"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "journal/recorder.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

// The benchmark is only run if --gtest_filter=PlayerTest.Benchmarks
void BM_PlayForReal(benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    Player player(
        R"(P:\Public Mockingbird\Principia\Journals\JOURNAL.20160626-143407)");
    int count = 0;
    while (player.Play()) {
      ++count;
      LOG_IF(ERROR, (count % 100'000) == 0)
          << count << " journal entries replayed";
    }
  }
}

BENCHMARK(BM_PlayForReal);

class PlayerTest : public ::testing::Test {
 protected:
  PlayerTest()
      : test_info_(testing::UnitTest::GetInstance()->current_test_info()),
        test_case_name_(test_info_->test_case_name()),
        test_name_(test_info_->name()),
        plugin_(interface::principia__NewPlugin("0 s", "0 s", 0)) {}

  template<typename Profile>
  bool RunIfAppropriate(serialization::Method const& method_in,
                        serialization::Method const& method_out_return,
                        Player& player) {
    return player.RunIfAppropriate<Profile>(method_in, method_out_return);
  }

  ::testing::TestInfo const* const test_info_;
  std::string const test_case_name_;
  std::string const test_name_;
  std::unique_ptr<ksp_plugin::Plugin> plugin_;
};

TEST_F(PlayerTest, PlayTiny) {
  // Do the recording in a separate thread to make sure that it activates using
  // a different static variable than the one in the plugin dynamic library.
  std::thread recorder([this]() {
    Recorder* const r(new Recorder(test_name_ + ".journal.hex"));
    Recorder::Activate(r);

    {
      Method<NewPlugin> m({"1 s", "2 s", 3});
      m.Return(plugin_.get());
    }
    {
      const ksp_plugin::Plugin* plugin = plugin_.get();
      Method<DeletePlugin> m({&plugin}, {&plugin});
      m.Return();
    }
    Recorder::Deactivate();
  });
  recorder.join();

  Player player(test_name_ + ".journal.hex");

  // Replay the journal.
  int count = 0;
  while (player.Play()) {
    ++count;
  }
  EXPECT_EQ(2, count);
}

// This test (a.k.a. benchmark) is only run if the --gtest_filter flag names it
// explicitly.
TEST_F(PlayerTest, Benchmarks) {
  if (testing::FLAGS_gtest_filter == test_case_name_ + "." + test_name_) {
    benchmark::RunSpecifiedBenchmarks();
  }
}

#if 1
// This test is only run if the --gtest_filter flag names it explicitly.
TEST_F(PlayerTest, Debug) {
  if (testing::FLAGS_gtest_filter == test_case_name_ + "." + test_name_) {
    // An example of how journaling may be used for debugging.  You must set
    // |path| and fill the |method_in| and |method_out_return| protocol buffers.
    std::string path =
        R"(P:\Public Mockingbird\Principia\Journals\JOURNAL.20170312-125704)";
    Player player(path);
    int count = 0;
    while (player.Play()) {
      ++count;
      LOG_IF(ERROR, (count % 100'000) == 0) << count
                                            << " journal entries replayed";
    }
    LOG(ERROR) << count << " journal entries in total";
    LOG(ERROR) << "Last successful method in:\n"
               << player.last_method_in().DebugString();
    LOG(ERROR) << "Last successful method out/return: \n"
               << player.last_method_out_return().DebugString();

#if 0
    serialization::Method method_in;
    auto* extension = method_in.MutableExtension(
        serialization::ReportCollision::extension);
    auto* in = extension->mutable_in();
    in->set_plugin(355375312);
    in->set_vessel1_guid("14b05bd3-9707-4d49-a6be-a7de481f3e0a");
    in->set_vessel2_guid("3e6fcb7e-4761-48ed-829f-0adb035f457e");
    serialization::Method method_out_return;
    method_out_return.MutableExtension(
        serialization::ReportCollision::extension);
    LOG(ERROR) << "Running unpaired method:\n" << method_in.DebugString();
    CHECK(RunIfAppropriate<ReportCollision>(method_in,
                                            method_out_return,
                                            player));
#endif
  }
}
#endif

}  // namespace journal
}  // namespace principia
