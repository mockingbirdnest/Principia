
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
    Player player(
        R"(P:\Public Mockingbird\Principia\Journals\JOURNAL.20160502-200332)");
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
        plugin_(interface::principia__NewPlugin(1, 2)),
        recorder_(new Recorder(test_name_ + ".journal.hex")) {
    Recorder::Activate(recorder_);
  }

  ~PlayerTest() override {
    Recorder::Deactivate();
  }

  template<typename Profile>
  void RunIfAppropriate(serialization::Method const& method_in,
                        serialization::Method const& method_out_return,
                        Player& player) {
    player.RunIfAppropriate<Profile>(method_in, method_out_return);
  }

  ::testing::TestInfo const* const test_info_;
  std::string const test_case_name_;
  std::string const test_name_;
  std::unique_ptr<ksp_plugin::Plugin> plugin_;
  Recorder* recorder_;
};

TEST_F(PlayerTest, PlayTiny) {
  {
    Method<NewPlugin> method_in({1, 2});
    method_in.Return(plugin_.get());
  }
  {
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> method_in({&plugin}, {&plugin});
    method_in.Return();
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

#if 0
// This test is only run if the --gtest_filter flag names it explicitly.
TEST_F(PlayerTest, Debug) {
  if (testing::FLAGS_gtest_filter == test_case_name_ + "." + test_name_) {
    // An example of how journalling may be used for debugging.  You must set
    // |path| and fill the |method_in| and |method_out_return| protocol buffers.
    std::string path =
        R"(P:\Public Mockingbird\Principia\Journals\JOURNAL.20160502-200332)";
    Player player(path);
    int count = 0;
    while (player.Play()) {
      ++count;
      LOG_IF(ERROR, (count % 100'000) == 0) << count
                                            << " journal entries replayed";
    }
    LOG(ERROR) << count << " journal entries in total";

    //serialization::Method method_in;
    //auto* extension = method_in.MutableExtension(
    //    serialization::SerializePlugin::extension);
    //auto* in = extension->mutable_in();
    //in->set_plugin(850673856);
    //in->set_serializer(0);
    //serialization::Method method_out_return;
    //method_out_return.MutableExtension(
    //    serialization::SerializePlugin::extension);
    //RunIfAppropriate<SerializePlugin>(method_in, method_out_return, player);
  }
}
#endif

}  // namespace journal
}  // namespace principia
