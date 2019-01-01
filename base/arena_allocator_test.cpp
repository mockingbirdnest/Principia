
#include "base/arena_allocator.hpp"

#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::Eq;
using ::testing::Ge;

class ArenaAllocatorTest : public ::testing::Test {
 protected:
  ArenaAllocatorTest() : allocator_(&arena_) {
    arena_.Init(google::protobuf::ArenaOptions());
  }

  google::protobuf::Arena arena_;
  ArenaAllocator<std::string> allocator_;
};

TEST_F(ArenaAllocatorTest, Container) {
  EXPECT_THAT(arena_.SpaceUsed(), Eq(0));
  std::vector<std::string, ArenaAllocator<std::string>> v(1000, allocator_);
  EXPECT_THAT(arena_.SpaceUsed(), Ge(1000 * sizeof(std::string)));
}

TEST_F(ArenaAllocatorTest, ContainerAndString) {
  EXPECT_THAT(arena_.SpaceUsed(), Eq(0));
  using S = std::basic_string<char,
                              std::char_traits<char>,
                              ArenaAllocator<char>>;
  ArenaAllocator<char> allocator1(&arena_);
  ArenaAllocator<S> allocator2(&arena_);
  std::vector<S, ArenaAllocator<S>> v(1000, S("foo", allocator1), allocator2);
  EXPECT_THAT(arena_.SpaceUsed(), Ge(1000 * (sizeof(std::string) + 4)));
}

}  // namespace base
}  // namespace principia
