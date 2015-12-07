#include "journal/player.hpp"

#include <list>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "journal/recorder.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

class PlayerTest : public testing::Test {
 protected:
  PlayerTest()
      : test_name_(
            testing::UnitTest::GetInstance()->current_test_info()->name()),
        plugin_(ksp_plugin::principia__NewPlugin(1, 2)),
        recorder_(new Recorder(test_name_ + ".journal.hex")) {
    Recorder::Activate(recorder_);
  }

  ~PlayerTest() override {
    Recorder::Deactivate();
  }

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

#if 0
TEST_F(PlayerTest, PlayForReal) {
  Player player("P:\\Public Mockingbird\\Principia\\JOURNAL.20151206-170008");
  int count = 0;
  while (player.Play()) {
    ++count;
    LOG_IF(ERROR, (count % 10'000) == 0)
        << count << " journal entries replayed";
  }
}
#endif

}  // namespace journal
}  // namespace principia
