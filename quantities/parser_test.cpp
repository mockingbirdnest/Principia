
#include "quantities/parser.hpp"

#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/astronomy.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using testing_utilities::AlmostEquals;

namespace quantities {

using astronomy::AstronomicalUnit;
using si::Day;
using si::Degree;
using si::Kilo;
using si::Metre;
using si::Radian;
using si::Second;
using si::Steradian;
using si::Watt;

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
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m ^ 3/s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m^3 / s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m ^ 3 / s ^ 2"));
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m^ 3/s ^2"));
}

TEST_F(ParserDeathTest, SpacesError) {
  EXPECT_DEATH({
    ParseQuantity<double>("1. 23");
  }, "Unsupported unit");
  EXPECT_DEATH({
    ParseQuantity<double>("1 .23");
  }, "Unsupported unit");
  EXPECT_DEATH({
    ParseQuantity<Speed>("1. 23 m/s");
  }, "Unsupported unit");
  EXPECT_DEATH({
    ParseQuantity<Speed>("1 .23 m/s");
  }, "Unsupported unit");
  EXPECT_DEATH({
    ParseQuantity<Speed>("1.23 m^- 2/s");
  }, "invalid integer");
}

TEST_F(ParserDeathTest, UnitError) {
  EXPECT_DEATH({
    ParseQuantity<Length>("1.23 nm");
  }, "Unsupported unit");
  EXPECT_DEATH({
    ParseQuantity<Time>("1.23 hr");
  }, "Unsupported unit");
  EXPECT_DEATH({
    ParseQuantity<Angle>("1.23 grd");
  }, "Unsupported unit");
  EXPECT_DEATH({
    ParseQuantity<Radiance>("1.23 W/sr m^2");
  }, "Unsupported unit sr m");
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
  EXPECT_EQ(1.23 * Degree, ParseQuantity<Angle>(u8"1.23 °"));
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
              AlmostEquals(1.23 * Pow<3>(AstronomicalUnit) / Pow<2>(Day), 2));
}

TEST_F(ParserTest, ParseAcceleration) {
  EXPECT_EQ(1.23 * Metre / Second / Second,
            ParseQuantity<Acceleration>("1.23 m/s/s"));
  EXPECT_EQ(1.23 * Kilo(Metre) / Second / Second,
            ParseQuantity<Acceleration>("1.23 km/s^2"));
}

TEST_F(ParserTest, ParseAngularFrequency) {
  EXPECT_EQ(1.23 * Radian / Second,
            ParseQuantity<AngularFrequency>("1.23 rad/s"));
  EXPECT_EQ(1.23 * Radian / Second,
            ParseQuantity<AngularFrequency>("1.23 s^ -1 rad"));
}

TEST_F(ParserTest, ParseRadiance) {
  EXPECT_EQ(1.23 * Watt / (Steradian * Metre * Metre),
    ParseQuantity<Radiance>("1.23 W sr^-1 m^-2"));
}

TEST_F(ParserTest, ParseFrequency) {
  EXPECT_EQ(1.23 / Second, ParseQuantity<Inverse<Time>>("1.23 / s"));
}

}  // namespace quantities
}  // namespace principia
