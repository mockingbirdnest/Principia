
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

void BM_PlayForReal(benchmark::State& state) {
  while (state.KeepRunning()) {
    Player player(
        R"(P:\Public Mockingbird\Principia\Journals\JOURNAL.20170422-143027)");
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

TEST_F(PlayerTest, DISABLED_Benchmarks) {
  benchmark::RunSpecifiedBenchmarks();
}

TEST_F(PlayerTest, DISABLED_Debug) {
  // An example of how journaling may be used for debugging.  You must set
  // |path| and fill the |method_in| and |method_out_return| protocol buffers.
  std::string path =
      R"(P:\Public Mockingbird\Principia\Journals\JOURNAL.20170422-143027)";
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

#if 1
  serialization::Method method_in;
  {
    auto* extension = method_in.MutableExtension(
        serialization::RenderedPredictionNodes::extension);
    auto* in = extension->mutable_in();
    in->set_plugin(3220795440);
    in->set_vessel_guid("cf6715cb-8780-4c39-8c01-3cc0cf35683b");
    in->mutable_sun_world_position()->set_x(-3412377724.9237709);
    in->mutable_sun_world_position()->set_y(5875636073.1931963);
    in->mutable_sun_world_position()->set_z(-11782002982.199373);
  }
  serialization::Method method_out_return;
  {
    auto* extension = method_out_return.MutableExtension(
        serialization::RenderedPredictionNodes::extension);
    auto* out = extension->mutable_out();
    out->set_ascending(1);
    out->set_descending(2);
  }
  LOG(ERROR) << "Running unpaired method:\n" << method_in.DebugString();
  CHECK(RunIfAppropriate<RenderedPredictionNodes>(method_in,
                                                  method_out_return,
                                                  player));
#endif
}

}  // namespace journal
}  // namespace principia
