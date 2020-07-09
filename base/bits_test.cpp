#include "base/bits.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace base {

TEST(BitsTest, FloorLog2) {
  EXPECT_EQ(0, FloorLog2(0));
  EXPECT_EQ(0, FloorLog2(1));
  EXPECT_EQ(1, FloorLog2(2));
  EXPECT_EQ(1, FloorLog2(3));
  EXPECT_EQ(2, FloorLog2(4));
  EXPECT_EQ(2, FloorLog2(7));
  EXPECT_EQ(3, FloorLog2(8));
}

TEST(BitsTest, PowerOf2Le) {
  EXPECT_EQ(0, PowerOf2Le(0));
  EXPECT_EQ(1, PowerOf2Le(1));
  EXPECT_EQ(2, PowerOf2Le(2));
  EXPECT_EQ(2, PowerOf2Le(3));
  EXPECT_EQ(4, PowerOf2Le(4));
  EXPECT_EQ(4, PowerOf2Le(7));
  EXPECT_EQ(8, PowerOf2Le(8));
}

}  // namespace base
}  // namespace principia
