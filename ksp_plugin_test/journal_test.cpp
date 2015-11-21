#include "ksp_plugin/journal.hpp"

#include <list>
#include <string>

#include "gtest/gtest.h"
#include "ksp_plugin/mock_plugin.hpp"
#include "serialization/journal.pb.h"

namespace principia {
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

TEST_F(JournalTest, Journal) {
  {
    const Plugin* plugin = plugin_.get();
    Journal::Method<DeletePlugin> m({plugin}, {&plugin});
    m.Return();
  }
  {
    Journal::Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
  }

  //auto const& j = *journal();
  //EXPECT_EQ(2, j.size());

  //auto it = j.begin();
  //{
  //  EXPECT_TRUE(it->HasExtension(serialization::DeletePlugin::extension));
  //  auto const& extension =
  //      it->GetExtension(serialization::DeletePlugin::extension);
  //  EXPECT_TRUE(extension.has_in());
  //  EXPECT_NE(0, extension.in().plugin());
  //  EXPECT_TRUE(extension.has_out());
  //  EXPECT_NE(0, extension.out().plugin());
  //  EXPECT_EQ(extension.in().plugin(), extension.out().plugin());
  //}
  //++it;
  //{
  //  EXPECT_TRUE(it->HasExtension(serialization::NewPlugin::extension));
  //  auto const& extension =
  //      it->GetExtension(serialization::NewPlugin::extension);
  //  EXPECT_TRUE(extension.has_in());
  //  EXPECT_EQ(1, extension.in().initial_time());
  //  EXPECT_EQ(2, extension.in().planetarium_rotation_in_degrees());
  //  EXPECT_TRUE(extension.has_return_());
  //  EXPECT_NE(0, extension.return_().plugin());
  //}
}

}  // namespace ksp_plugin
}  // namespace principia
