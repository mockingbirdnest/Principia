#include "quantities/parser.hpp"

#include "gtest/gtest.h"

#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using si::AstronomicalUnit;
using si::Kilo;
using si::Metre;
using si::Second;

namespace quantities {

class ParserTest : public ::testing::Test {
};

using ParserDeathTest = ParserTest;

TEST_F(ParserTest, SpacesSuccess) {
  EXPECT_EQ(1.23, ParseQuantity<double>("  1.23"));
  EXPECT_EQ(1.23, ParseQuantity<double>("  1.23   "));
  EXPECT_EQ(1.23, ParseQuantity<double>("1.23 "));
  EXPECT_EQ(1.23, ParseQuantity<double>("1.23"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("  1.23 m/s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("  1.23 m/s   "));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/s "));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m /s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/ s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m / s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23m/s"));
}

TEST_F(ParserDeathTest, SpacesError) {
  EXPECT_DEATH({
    ParseQuantity<double>("1. 23");
  }, "empty");
  EXPECT_DEATH({
    ParseQuantity<double>("1 .23");
  }, "empty");
  EXPECT_DEATH({
    ParseQuantity<Speed>("1. 23 m/s");
  }, "Unsupported.*length");
  EXPECT_DEATH({
    ParseQuantity<Speed>("1 .23 m/s");
  }, "Unsupported.*length");
}

TEST_F(ParserDeathTest, UnitError) {
  EXPECT_DEATH({
    ParseQuantity<Length>("1.23 nm");
  }, "Unsupported.*length");
}

TEST_F(ParserTest, ParseDouble) {
  EXPECT_EQ(1.23, ParseQuantity<double>("1.23"));
  EXPECT_EQ(-3.45, ParseQuantity<double>("-3.45"));
}

TEST_F(ParserTest, ParseLength) {
  EXPECT_EQ(1.23 * Metre, ParseQuantity<Length>("1.23 m"));
  EXPECT_EQ(1.23 * Kilo(Metre), ParseQuantity<Length>("1.23 km"));
  EXPECT_EQ(1.23 * AstronomicalUnit, ParseQuantity<Length>("1.23 au"));
}

}  // namespace quantities
}  // namespace principia
