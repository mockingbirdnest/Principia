
#include "journal/player.hpp"

#include <list>
#include <string>
#include <vector>

#include "benchmark/benchmark.h"
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
    Player player(R"(P:\Public Mockingbird\Principia\JOURNAL.20160207-192649)");
    int count = 0;
    while (player.Play()) {
      ++count;
      LOG_IF(ERROR, (count % 100'000) == 0)
          << count << " journal entries replayed";
    }
  }
}

BENCHMARK(BM_PlayForReal);

#if 0
// An example of how journalling may be used for debugging.  You must set |path|
// and fill the |m| protocol buffer.  This probably requires to make
// RunIfAppropriate public.
void BM_PlayForDebug(benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    std::string path;
    Player player(path);
    int count = 0;
    while (player.Play()) {
      ++count;
      LOG_IF(ERROR, (count % 100'000) == 0)
          << count << " journal entries replayed";
    }

    serialization::Method m;
    auto* extension =
        m.MutableExtension(serialization::FlightPlanGetManoeuvre::extension);
    auto* in = extension->mutable_in();
    in->set_plugin(1093188824);
    in->set_vessel_guid("32100de3-7be2-432e-ac44-d34d2890bcea");
    in->set_index(0);
    player.RunIfAppropriate<FlightPlanGetManoeuvre>(m);
  }
}

BENCHMARK(BM_PlayForDebug);
#endif

class PlayerTest : public ::testing::Test {
 protected:
  PlayerTest()
      : test_info_(testing::UnitTest::GetInstance()->current_test_info()),
        test_case_name_(test_info_->test_case_name()),
        test_name_(test_info_->name()),
        plugin_(interface::principia__NewPlugin(1, 2)),
        recorder_(new Recorder(test_name_ + ".journal.hex",
                               /*verbose=*/false)) {
    Recorder::Activate(recorder_);
  }

  ~PlayerTest() override {
    Recorder::Deactivate();
  }

  ::testing::TestInfo const* const test_info_;
  std::string const test_case_name_;
  std::string const test_name_;
  std::unique_ptr<ksp_plugin::Plugin> plugin_;
  Recorder* recorder_;
};

TEST_F(PlayerTest, PlayTiny) {
  {
    Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
  }
  {
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> m({&plugin}, {&plugin});
    m.Return();
  }

  Player player(test_name_ + ".journal.hex");

  // Replay the journal.  Note that the journal doesn't grow as we replay
  // because we didn't call principia__ActivateRecorder so there is no active
  // journal in the ksp_plugin assembly.
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

}  // namespace journal
}  // namespace principia
