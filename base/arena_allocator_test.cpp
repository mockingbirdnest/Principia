
#include "base/arena_allocator.hpp"

#include <string>

#include "gtest/gtest.h"

namespace principia {
namespace base {

class ArenaAllocatorTest : public ::testing::Test {
 protected:
  ArenaAllocatorTest() : allocator_(&arena_) {
    arena_.Init(google::protobuf::ArenaOptions());
  }

  google::protobuf::Arena arena_;
  ArenaAllocator<std::string> allocator_;
};

TEST_F(ArenaAllocatorTest, Container) {
  CHECK_EQ(0, arena_.SpaceUsed());
  std::vector<std::string, ArenaAllocator<std::string>> v(1000, allocator_);
  CHECK_LT(0, arena_.SpaceUsed());
}

TEST_F(ArenaAllocatorTest, PlacementNew) {
  CHECK_EQ(0, arena_.SpaceUsed());
  std::string* s = new(allocator_) std::string("foo");
  CHECK_LT(0, arena_.SpaceUsed());
}

}  // namespace base
}  // namespace principia
