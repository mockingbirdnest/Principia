#include "base/flags.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace base {

TEST(FlagsTest, Basics) {
  Flags::Set("zfp=yes,hex,gipfeli=");
  EXPECT_TRUE(Flags::IsPresent("zfp"));
  EXPECT_TRUE(Flags::IsPresent("hex"));
  EXPECT_TRUE(Flags::IsPresent("gipfeli"));
  EXPECT_FALSE(Flags::IsPresent("decimal"));
  EXPECT_EQ("yes", Flags::Value("zfp"));
  EXPECT_EQ("", Flags::Value("hex"));
  EXPECT_EQ("", Flags::Value("gipfeli"));
}

TEST(FlagsTest, Space) {
  Flags::Set("zfp=yes , hex");
  EXPECT_TRUE(Flags::IsPresent("zfp"));
  EXPECT_TRUE(Flags::IsPresent(" hex"));
  EXPECT_EQ("yes", Flags::Value("zfp"));
  EXPECT_EQ("", Flags::Value(" hex"));
}

}  // namespace base
}  // namespace principia
