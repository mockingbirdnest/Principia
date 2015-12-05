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
        recorder_(new Recorder(test_name_)) {
    Recorder::Activate(recorder_);
  }

  ~PlayerTest() override {
    Recorder::Deactivate();
  }

  static std::vector<serialization::Method> ReadAll(
      std::experimental::filesystem::path const& path) {
    std::vector<serialization::Method> methods;
    Player player(path);
    for (std::unique_ptr<serialization::Method> method = player.Read();
         method != nullptr;
         method = player.Read()) {
      methods.push_back(*method);
    }
    return methods;
  }


  std::string const test_name_;
  std::unique_ptr<ksp_plugin::Plugin> plugin_;
  Recorder* recorder_;
};

TEST_F(PlayerTest, Playing) {
  {
    Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
  }
  {
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> m({&plugin}, {&plugin});
    m.Return();
  }

  // Read all the messages to determine the current size of the journal.
  std::vector<serialization::Method> methods1 = ReadAll(test_name_);
  Player player(test_name_);

  // Replay the journal.  Note that the journal doesn't grow as we replay
  // because we didn't call principia__ActivateRecorder so there is no active
  // journal in the ksp_plugin assembly.
  while (player.Play()) {}
}

}  // namespace journal
}  // namespace principia
