#include "journal/player.hpp"

#include <chrono>
#include <memory>
#include <string>
#include <thread>

#include "benchmark/benchmark.h"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "journal/recorder.hpp"
#include "ksp_plugin/plugin.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

using namespace principia::journal::_method;
using namespace principia::journal::_player;
using namespace principia::journal::_recorder;
using namespace principia::ksp_plugin::_plugin;
using namespace std::chrono_literals;

void BM_PlayForReal(benchmark::State& state) {
  for (auto _ : state) {
    Player player(
        R"(P:\Public Mockingbird\Principia\Journals\JOURNAL.20180311-192733)");
    int count = 0;
    while (player.Play(count)) {
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
        test_name_(test_info_->name()) {}

  template<typename Profile>
  bool RunIfAppropriate(serialization::Method const& method_in,
                        serialization::Method const& method_out_return,
                        Player& player) {
    return player.RunIfAppropriate<Profile>(method_in, method_out_return);
  }

  ::testing::TestInfo const* const test_info_;
  std::string const test_case_name_;
  std::string const test_name_;
};

TEST_F(PlayerTest, PlayTiny) {
  {
    std::unique_ptr<Plugin> const plugin(
        interface::principia__NewPlugin("MJD0", "MJD0", 0));

    Recorder* const r(new Recorder(test_name_ + ".journal.hex"));
    Recorder::Activate(r);

    {
      Method<NewPlugin> m({"MJD1", "MJD2", 3});
      m.Return(plugin.get());
    }
    {
      const Plugin* p = plugin.get();
      Method<DeletePlugin> m({&p}, {&p});
      m.Return();
    }
    Recorder::Deactivate();
  }

  Player player(test_name_ + ".journal.hex");

  // Replay the journal.
  int count = 0;
  while (player.Play(count)) {
    ++count;
  }
  EXPECT_EQ(3, count);
}

TEST_F(PlayerTest, DISABLED_SECULAR_Benchmarks) {
  benchmark::RunSpecifiedBenchmarks();
}

// A convenience test to find the last unpaired method of a journal.  You must
// set `path`.
TEST_F(PlayerTest, DISABLED_SECULAR_Scan) {
  std::string path =
      R"(P:\Public Mockingbird\Principia\Crashes\3375\JOURNAL.20220610-092143)";  // NOLINT
  Player player(path);
  int count = 0;
  while (player.Scan(count)) {
    ++count;
    LOG_IF(ERROR, (count % 100'000) == 0) << count
                                          << " journal entries replayed";
  }
  LOG(ERROR) << count << " journal entries in total";
  LOG(ERROR) << "Last successful method in:\n"
             << player.last_method_in().DebugString();
  LOG(ERROR) << "Last successful method out/return: \n"
             << player.last_method_out_return().DebugString();
}

// A test to debug a journal.  You must set `path` and fill the `method_in` and
// `method_out_return` protocol buffers.
TEST_F(PlayerTest, DISABLED_SECULAR_Debug) {
  std::string path =
      R"(P:\Public Mockingbird\Principia\Issues\3872\JOURNAL.20240210-173425)";  // NOLINT
  Player player(path);
  int count = 0;
  while (player.Play(count)) {
    ++count;
    // Reset logging after each method so as to output all messages irrespective
    // of what the game did.
    google::LogToStderr();
    LOG_IF(ERROR, (count % 100'000) == 0) << count
                                          << " journal entries replayed";
  }
  LOG(ERROR) << count << " journal entries in total";
  LOG(ERROR) << "Last successful method in:\n"
             << player.last_method_in().DebugString();
  LOG(ERROR) << "Last successful method out/return: \n"
             << player.last_method_out_return().DebugString();
  std::this_thread::sleep_for(10s);

#if 0
  serialization::Method method_in;
  {
    auto* extension = method_in.MutableExtension(
        serialization::CollisionDeleteExecutor::extension);
    auto* in = extension->mutable_in();
    in->set_plugin(2237555212240);
    in->set_executor(2237561081696);
  }
  serialization::Method method_out_return;
  {
    auto* extension = method_out_return.MutableExtension(
        serialization::CollisionDeleteExecutor::extension);
  }
  LOG(ERROR) << "Running unpaired method:\n" << method_in.DebugString();
  CHECK(RunIfAppropriate<CollisionDeleteExecutor>(
      method_in, method_out_return, player));
#endif
#if 0
  std::this_thread::sleep_for(10s);
#endif
}

}  // namespace journal
}  // namespace principia
