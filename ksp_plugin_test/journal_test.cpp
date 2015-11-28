#include "ksp_plugin/journal.hpp"

#include <list>
#include <string>
#include <vector>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "gtest/gtest.h"
#include "ksp_plugin/mock_plugin.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::HexadecimalDecode;
using base::UniqueBytes;

namespace ksp_plugin {

class JournalTest : public testing::Test {
 protected:
  JournalTest()
      : test_name_(
            testing::UnitTest::GetInstance()->current_test_info()->name()),
        plugin_(std::make_unique<MockPlugin>()),
        journal_(new Journal(test_name_)) {
    Journal::Activate(journal_);
  }

  ~JournalTest() override {
    Journal::Deactivate();
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
  std::unique_ptr<MockPlugin> plugin_;
  Journal* journal_;
};

using JournalDeathTest = JournalTest;

TEST_F(JournalDeathTest, Return) {
  EXPECT_DEATH({
    Journal::Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
    m.Return();
  },
  "!returned_");
  EXPECT_DEATH({
    const Plugin* plugin = plugin_.get();
    Journal::Method<DeletePlugin> m({plugin}, {&plugin});
    m.Return();
    m.Return();
  },
  "!returned_");
  EXPECT_DEATH({
    Journal::Method<NewPlugin> m({1, 2});
  },
  "returned_");
}

TEST_F(JournalTest, Recording) {
  {
    const Plugin* plugin = plugin_.get();
    Journal::Method<DeletePlugin> m({plugin}, {&plugin});
    m.Return();
  }
  {
    Journal::Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
  }

  std::vector<serialization::Method> methods = ReadAll(test_name_);
  EXPECT_EQ(2, methods.size());
  auto it = methods.begin();
  {
    EXPECT_TRUE(it->HasExtension(serialization::DeletePlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::DeletePlugin::extension);
    EXPECT_TRUE(extension.has_in());
    EXPECT_NE(0, extension.in().plugin());
    EXPECT_TRUE(extension.has_out());
    EXPECT_NE(0, extension.out().plugin());
    EXPECT_EQ(extension.in().plugin(), extension.out().plugin());
  }
  ++it;
  {
    EXPECT_TRUE(it->HasExtension(serialization::NewPlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::NewPlugin::extension);
    EXPECT_TRUE(extension.has_in());
    EXPECT_EQ(1, extension.in().initial_time());
    EXPECT_EQ(2, extension.in().planetarium_rotation_in_degrees());
    EXPECT_TRUE(extension.has_return_());
    EXPECT_NE(0, extension.return_().new_plugin());
  }
}

TEST_F(JournalTest, Playing) {
  {
    Journal::Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
  }
  {
    const Plugin* plugin = plugin_.get();
    Journal::Method<DeletePlugin> m({plugin}, {&plugin});
    m.Return();
  }

  // Read all the messages to determine the current size of the journal.
  std::vector<serialization::Method> methods1 = ReadAll(test_name_);
  Player player(test_name_);

  // Replay the part of the journal written so far.  Cannot use the boolean
  // returned by Play because the journal grows as we replay.
  for (int i = 0; i < methods1.size(); ++i) {
    player.Play();
  }

  // Check that the journal has grown as expected.
  std::vector<serialization::Method> methods2 = ReadAll(test_name_);
  CHECK_EQ(2 * methods1.size(), methods2.size());
}

}  // namespace ksp_plugin
}  // namespace principia
