
#include "journal/player.hpp"

#include <list>
#include <string>
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

void BM_PlayForReal(benchmark::State& state) {
  while (state.KeepRunning()) {
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
        test_name_(test_info_->name()),
        plugin_(interface::principia__NewPlugin("MJD0", "MJD0", 0)) {}

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
  {
    Recorder* const r(new Recorder(test_name_ + ".journal.hex"));
    Recorder::Activate(r);

    {
      Method<NewPlugin> m({"MJD1", "MJD2", 3});
      m.Return(plugin_.get());
    }
    {
      const ksp_plugin::Plugin* plugin = plugin_.get();
      Method<DeletePlugin> m({&plugin}, {&plugin});
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

TEST_F(PlayerTest, DISABLED_SECULAR_Debug) {
  google::LogToStderr();
  // An example of how journaling may be used for debugging.  You must set
  // |path| and fill the |method_in| and |method_out_return| protocol buffers.
  std::string path =
      R"(P:\Public Mockingbird\Principia\Crashes\2507\JOURNAL.20200328-214524)";  // NOLINT
  Player player(path);
  int count = 0;
  while (player.Play(count)) {
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
  {
    auto* extension = method_in.MutableExtension(
        serialization::CatchUpLaggingVessels::extension);
    auto* in = extension->mutable_in();
    in->set_plugin(2014156823824);
  }
  serialization::Method method_out_return;
  {
    auto* extension = method_out_return.MutableExtension(
        serialization::CatchUpLaggingVessels::extension);
  }
  LOG(ERROR) << "Running unpaired method:\n" << method_in.DebugString();
  CHECK(RunIfAppropriate<CatchUpLaggingVessels>(
      method_in, method_out_return, player));
#endif
}

}  // namespace journal
}  // namespace principia
