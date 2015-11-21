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

  // TODO(phl): Some of this code should be moved to the methods that replay a
  // journal (TBD).
  std::vector<serialization::Method> methods;
  std::ifstream stream(test_name_, std::ios::in);
  std::string line(1000, ' ');
  while (stream.getline(&line[0], line.size())) {
    uint8_t const* const hexadecimal =
        reinterpret_cast<uint8_t const*>(line.c_str());
    int const hexadecimal_size = strlen(line.c_str());
    UniqueBytes bytes(hexadecimal_size >> 1);
    HexadecimalDecode({hexadecimal, hexadecimal_size},
                      {bytes.data.get(), bytes.size});
    methods.emplace_back();
    CHECK(methods.back().ParseFromArray(bytes.data.get(),
                                        static_cast<int>(bytes.size)));
  }

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
    EXPECT_NE(0, extension.return_().plugin());
  }
}

}  // namespace ksp_plugin
}  // namespace principia
