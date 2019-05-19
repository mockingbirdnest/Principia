
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
      R"(P:\Public Mockingbird\Principia\Crashes\2173\JOURNAL.20190519-024654)";  // NOLINT
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
        serialization::SetPartApparentDegreesOfFreedom::extension);
    auto* in = extension->mutable_in();
    in->set_plugin(6351018528);
    in->set_part_id(0);
    {
      auto* const dof = in->mutable_degrees_of_freedom();
      auto* const q = dof->mutable_q();
      q->set_x(-0.36689972877502441);
      q->set_y(-0.21265912055969238);
      q->set_z(-0.13075445592403412);
      auto* const p = dof->mutable_p();
      p->set_x(0.21613305807113647);
      p->set_y(0.12690158188343048);
      p->set_z(0.08301004022359848);
    }
    {
      auto* const dof = in->mutable_main_body_degrees_of_freedom();
      auto* const q = dof->mutable_q();
      q->set_x(-5269720.6347961426);
      q->set_y(-3050611.0277061462);
      q->set_z(-1874750.4663391113);
      auto* const p = dof->mutable_p();
      p->set_x(-0.0);
      p->set_y(-0.0);
      p->set_z(-0.0);
    }
  }
  serialization::Method method_out_return;
  {
    auto* extension = method_out_return.MutableExtension(
        serialization::SetPartApparentDegreesOfFreedom::extension);
  }
  LOG(ERROR) << "Running unpaired method:\n" << method_in.DebugString();
  CHECK(RunIfAppropriate<SetPartApparentDegreesOfFreedom>(method_in,
                                                          method_out_return,
                                                          player));
#endif
}

}  // namespace journal
}  // namespace principia
