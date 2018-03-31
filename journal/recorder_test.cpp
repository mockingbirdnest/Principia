
#include "journal/recorder.hpp"

#include <filesystem>
#include <list>
#include <string>
#include <vector>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "gtest/gtest.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

class RecorderTest : public testing::Test {
 protected:
  RecorderTest()
      : test_name_(
            testing::UnitTest::GetInstance()->current_test_info()->name()),
        plugin_(interface::principia__NewPlugin("MJD0", "MJD0", 0)),
        recorder_(new Recorder(test_name_ + ".journal.hex")) {
    Recorder::Activate(recorder_);
  }

  ~RecorderTest() override {
    Recorder::Deactivate();
  }

  static std::vector<serialization::Method> ReadAll(
      std::filesystem::path const& path) {
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

using JournalDeathTest = RecorderTest;

TEST_F(JournalDeathTest, Return) {
  EXPECT_DEATH({
    Method<NewPlugin> m({"1 s", "2 s", 3});
    m.Return(plugin_.get());
    m.Return(plugin_.get());
  },
  "!returned_");
  EXPECT_DEATH({
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> m({&plugin}, {&plugin});
    m.Return();
    m.Return();
  },
  "!returned_");
  EXPECT_DEATH({
    Method<NewPlugin> m({"1 s", "2 s", 3});
  },
  "returned_");
}

TEST_F(RecorderTest, Recording) {
  {
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> m({&plugin}, {&plugin});
    m.Return();
  }
  {
    Method<NewPlugin> m({"1 s", "2 s", 3});
    m.Return(plugin_.get());
  }

  std::vector<serialization::Method> const methods =
      ReadAll(test_name_ + ".journal.hex");
  EXPECT_EQ(4, methods.size());
  auto it = methods.begin();
  {
    EXPECT_TRUE(it->HasExtension(serialization::DeletePlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::DeletePlugin::extension);
    EXPECT_TRUE(extension.has_in());
    EXPECT_NE(0, extension.in().plugin());
    EXPECT_FALSE(extension.has_out());
  }
  ++it;
  {
    EXPECT_TRUE(it->HasExtension(serialization::DeletePlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::DeletePlugin::extension);
    EXPECT_FALSE(extension.has_in());
    EXPECT_TRUE(extension.has_out());
    EXPECT_NE(0, extension.out().plugin());
  }
  ++it;
  {
    EXPECT_TRUE(it->HasExtension(serialization::NewPlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::NewPlugin::extension);
    EXPECT_TRUE(extension.has_in());
    EXPECT_EQ("1 s", extension.in().game_epoch());
    EXPECT_EQ("2 s", extension.in().solar_system_epoch());
    EXPECT_EQ(3, extension.in().planetarium_rotation_in_degrees());
    EXPECT_FALSE(extension.has_return_());
  }
  ++it;
  {
    EXPECT_TRUE(it->HasExtension(serialization::NewPlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::NewPlugin::extension);
    EXPECT_FALSE(extension.has_in());
    EXPECT_TRUE(extension.has_return_());
    EXPECT_NE(0, extension.return_().result());
  }
}

}  // namespace journal
}  // namespace principia
