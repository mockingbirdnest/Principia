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

TEST(BitsTest, BitReversedIncrement) {
  EXPECT_EQ(0x8, BitReversedIncrement(0x0, 4));
  EXPECT_EQ(0x9, BitReversedIncrement(0x1, 4));
  EXPECT_EQ(0xA, BitReversedIncrement(0x2, 4));
  EXPECT_EQ(0xB, BitReversedIncrement(0x3, 4));
  EXPECT_EQ(0xC, BitReversedIncrement(0x4, 4));
  EXPECT_EQ(0xD, BitReversedIncrement(0x5, 4));
  EXPECT_EQ(0xE, BitReversedIncrement(0x6, 4));
  EXPECT_EQ(0xF, BitReversedIncrement(0x7, 4));
  EXPECT_EQ(0x4, BitReversedIncrement(0x8, 4));
  EXPECT_EQ(0x5, BitReversedIncrement(0x9, 4));
  EXPECT_EQ(0x6, BitReversedIncrement(0xA, 4));
  EXPECT_EQ(0x7, BitReversedIncrement(0xB, 4));
  EXPECT_EQ(0x2, BitReversedIncrement(0xC, 4));
  EXPECT_EQ(0x3, BitReversedIncrement(0xD, 4));
  EXPECT_EQ(0x1, BitReversedIncrement(0xE, 4));

}

}  // namespace base
}  // namespace principia
