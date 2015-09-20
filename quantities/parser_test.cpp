#include "quantities/parser.hpp"

#include "gtest/gtest.h"

#include "quantities/si.hpp"

namespace principia {

using si::Metre;

namespace quantities {

class ParserTest : public ::testing::Test {
};

TEST_F(ParserTest, ParserLength) {
  EXPECT_EQ(1.23 * Metre, ParseQuantity<Length>(" 1.23 m "));
}

}  // namespace quantities
}  // namespace principia
