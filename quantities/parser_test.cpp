#include "quantities/parser.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using si::AstronomicalUnit;
using si::Day;
using si::Degree;
using si::Kilo;
using si::Metre;
using si::Radian;
using si::Second;
using testing_utilities::AlmostEquals;

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

TEST_F(ParserTest, ParseAngle) {
  EXPECT_EQ(1.23 * Degree, ParseQuantity<Angle>("1.23 deg"));
  EXPECT_EQ(1.23 * Degree, ParseQuantity<Angle>("1.23 °"));
  EXPECT_EQ(1.23 * Radian, ParseQuantity<Angle>("1.23 rad"));
}

TEST_F(ParserTest, ParseSpeed) {
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/s"));
  EXPECT_EQ(1.23 * Kilo(Metre) / Second, ParseQuantity<Speed>("1.23 km/s"));
  EXPECT_EQ(1.23 * Kilo(Metre) / Day, ParseQuantity<Speed>("1.23 km/d"));
  EXPECT_EQ(1.23 * AstronomicalUnit / Day, ParseQuantity<Speed>("1.23 au/d"));
}

TEST_F(ParserTest, ParseGravitationalParameter) {
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m^3/s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Kilo(Metre)) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 km^3/s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Kilo(Metre)) / Pow<2>(Day),
            ParseQuantity<GravitationalParameter>("1.23 km^3/d^2"));
  EXPECT_THAT(ParseQuantity<GravitationalParameter>("1.23 au^3/d^2"),
              AlmostEquals(1.23 * Pow<3>(AstronomicalUnit) / Pow<2>(Day), 1));
}

}  // namespace quantities
}  // namespace principia
