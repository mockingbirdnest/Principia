#include "base/flags.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::ElementsAre;
using ::testing::IsEmpty;

TEST(FlagsTest, Basics) {
  Flags::Set("zfp", "yes");
  Flags::Set("hex", "");
  Flags::Set("gipfeli", "on");
  Flags::Set("gipfeli", "off");
  EXPECT_TRUE(Flags::IsPresent("zfp"));
  EXPECT_TRUE(Flags::IsPresent("hex"));
  EXPECT_TRUE(Flags::IsPresent("gipfeli"));
  EXPECT_FALSE(Flags::IsPresent("decimal"));
  EXPECT_TRUE(Flags::IsPresent("zfp", "yes"));
  EXPECT_TRUE(Flags::IsPresent("hex", ""));
  EXPECT_TRUE(Flags::IsPresent("gipfeli", "on"));
  EXPECT_TRUE(Flags::IsPresent("gipfeli", "off"));
  EXPECT_FALSE(Flags::IsPresent("gipfeli", "n/a"));
  EXPECT_FALSE(Flags::IsPresent("decimal", "0"));
  EXPECT_THAT(Flags::Values("zfp"), ElementsAre("yes"));
  EXPECT_THAT(Flags::Values("hex"), ElementsAre(""));
  EXPECT_THAT(Flags::Values("gipfeli"), ElementsAre("off", "on"));
  EXPECT_THAT(Flags::Values("decimal"), IsEmpty());
}

}  // namespace base
}  // namespace principia
