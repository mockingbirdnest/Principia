
#include <functional>
#include <string>

#include "google/protobuf/stubs/common.h"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/astronomy.hpp"
#include "quantities/bipm.hpp"
#include "quantities/constants.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/algebra.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace quantities {

using astronomy::EarthMass;
using astronomy::JulianYear;
using astronomy::JupiterMass;
using astronomy::LightYear;
using astronomy::LunarDistance;
using astronomy::Parsec;
using astronomy::SolarMass;
using constants::ElectronMass;
using constants::GravitationalConstant;
using constants::SpeedOfLight;
using constants::StandardGravity;
using constants::VacuumPermeability;
using constants::VacuumPermittivity;
using si::Ampere;
using si::AstronomicalUnit;
using si::Candela;
using si::Cycle;
using si::Day;
using si::Degree;
using si::Hour;
using si::Kelvin;
using si::Kilogram;
using si::Mega;
using si::Metre;
using si::Mole;
using si::Radian;
using si::Second;
using si::Steradian;
using testing_utilities::AlmostEquals;
using testing_utilities::RelativeError;
using testing_utilities::VanishesBefore;
using uk::Foot;
using uk::Furlong;
using uk::Mile;
using uk::Rood;
using ::testing::Lt;
using ::testing::MatchesRegex;

class QuantitiesTest : public testing::Test {
 protected:
};

using QuantitiesDeathTest = QuantitiesTest;

TEST_F(QuantitiesTest, AbsoluteValue) {
  EXPECT_EQ(Abs(-1729), 1729);
  EXPECT_EQ(Abs(1729), 1729);
}

TEST_F(QuantitiesTest, DimensionfulComparisons) {
  testing_utilities::TestOrder(EarthMass, JupiterMass);
  testing_utilities::TestOrder(LightYear, Parsec);
  testing_utilities::TestOrder(-SpeedOfLight, SpeedOfLight);
  testing_utilities::TestOrder(SpeedOfLight * Day, LightYear);
}

TEST_F(QuantitiesTest, DimensionlfulOperations) {
  testing_utilities::TestVectorSpace(
      0 * Metre / Second, SpeedOfLight, 88 * Mile / Hour,
      -340.29 * Metre / Second, 0.0, 1.0, -2 * π, 1729.0, 0, 2);
  // Dimensionful multiplication is a tensor product, see [Tao 2012].
  testing_utilities::TestBilinearMap(
      std::multiplies<>(), SolarMass,
      ElectronMass, SpeedOfLight, 1 * Furlong / JulianYear, -e, 0, 2);
}

TEST_F(QuantitiesTest, DimensionlessExponentiation) {
  double const number   = π - 42;
  double positivePower = 1;
  double negativePower = 1;
  EXPECT_EQ(positivePower, Pow<0>(number));
  positivePower *= number;
  negativePower /= number;
  EXPECT_EQ(positivePower, Pow<1>(number));
  EXPECT_EQ(negativePower, Pow<-1>(number));
  positivePower *= number;
  negativePower /= number;
  EXPECT_EQ(positivePower, Pow<2>(number));
  EXPECT_EQ(negativePower, Pow<-2>(number));
  positivePower *= number;
  negativePower /= number;
  EXPECT_EQ(positivePower, Pow<3>(number));
  EXPECT_EQ(negativePower, Pow<-3>(number));
  positivePower *= number;
  negativePower /= number;
  // This one calls |std::pow|.
  EXPECT_THAT(positivePower, AlmostEquals(Pow<4>(number), 0, 1));
  EXPECT_THAT(negativePower, AlmostEquals(Pow<-4>(number), 0, 1));
}

// The Greek letters cause a warning when stringified by the macros, because
// apparently Visual Studio doesn't encode strings in UTF-8 by default.
#pragma warning(disable: 4566)

TEST_F(QuantitiesTest, Formatting) {
  auto const all_the_units = 1 * Metre * Kilogram * Second * Ampere * Kelvin /
                                 (Mole * Candela * Cycle * Radian * Steradian);
  std::string const expected = std::string("+1e+00 m kg s A K mol^-1") +
                               " cd^-1 cycle^-1 rad^-1 sr^-1";
  std::string const actual = DebugString(all_the_units, 0);
  EXPECT_EQ(expected, actual);
  std::string const π17 = "\\+3\\.1415926535897931.e\\+00";
  EXPECT_THAT(DebugString(π), MatchesRegex(π17));
  std::string const minus_e17 = "\\-2\\.718281828459045..e\\+00";
  EXPECT_THAT(DebugString(-e), MatchesRegex(minus_e17));
}

TEST_F(QuantitiesTest, PhysicalConstants) {
  // By definition.
  EXPECT_THAT(1 / Pow<2>(SpeedOfLight),
              AlmostEquals(VacuumPermittivity * VacuumPermeability, 1));
  // The Keplerian approximation for the mass of the Sun
  // is fairly accurate.
  EXPECT_THAT(RelativeError(
                  4 * Pow<2>(π) * Pow<3>(AstronomicalUnit) /
                      (GravitationalConstant * Pow<2>(JulianYear)),
                  SolarMass),
              Lt(4E-5));
  EXPECT_THAT(RelativeError(1 * Parsec, 3.26156 * LightYear), Lt(2E-6));
  // The Keplerian approximation for the mass of the Earth
  // is pretty bad, but the error is still only 1%.
  EXPECT_THAT(RelativeError(
                  4 * Pow<2>(π) * Pow<3>(LunarDistance) /
                      (GravitationalConstant * Pow<2>(27.321582 * Day)),
                  EarthMass),
              Lt(1E-2));
  EXPECT_THAT(RelativeError(1 * SolarMass, 1047 * JupiterMass), Lt(4E-4));
  // Delambre & Méchain.
  EXPECT_THAT(RelativeError(
                  GravitationalConstant * EarthMass /
                      Pow<2>(40 * Mega(Metre) / (2 * π)),
                  StandardGravity),
              Lt(4E-3));
  // Talleyrand.
  EXPECT_THAT(RelativeError(π * Sqrt(1 * Metre / StandardGravity), 1 * Second),
              Lt(4E-3));
}

#pragma warning(default: 4566)

TEST_F(QuantitiesTest, TrigonometricFunctions) {
  EXPECT_EQ(Cos(0 * Degree), 1);
  EXPECT_EQ(Sin(0 * Degree), 0);
  EXPECT_THAT(Cos(90 * Degree), VanishesBefore(1.0, 0));
  EXPECT_EQ(Sin(90 * Degree), 1);
  EXPECT_EQ(Cos(180 * Degree), -1);
  EXPECT_THAT(Sin(180 * Degree), VanishesBefore(1.0, 1));
  EXPECT_THAT(Cos(-90 * Degree), VanishesBefore(1.0, 0));
  EXPECT_EQ(Sin(-90 * Degree), -1);
  for (int k = 1; k < 360; ++k) {
    // Don't test for multiples of 90 degrees as zeros lead to horrible
    // conditioning.
    if (k % 90 != 0) {
      EXPECT_THAT(Cos((90 - k) * Degree),
                  AlmostEquals(Sin(k * Degree), 0, 47));
      EXPECT_THAT(Sin(k * Degree) / Cos(k * Degree),
                  AlmostEquals(Tan(k * Degree), 0, 2));
      EXPECT_THAT(((k + 179) % 360 - 179) * Degree,
                  AlmostEquals(ArcTan(Sin(k * Degree), Cos(k * Degree)),
                               0, 77));
      EXPECT_THAT(((k + 179) % 360 - 179) * Degree,
                  AlmostEquals(ArcTan(Sin(k * Degree) * AstronomicalUnit,
                                      Cos(k * Degree) * AstronomicalUnit),
                               0, 77));
      EXPECT_THAT(Cos(ArcCos(Cos(k * Degree))),
                  AlmostEquals(Cos(k * Degree), 0, 7));
      EXPECT_THAT(Sin(ArcSin(Sin(k * Degree))),
                  AlmostEquals(Sin(k * Degree), 0, 1));
    }
  }
  // Horribly conditioned near 0, so not in the loop above.
  EXPECT_EQ(Tan(ArcTan(Tan(-42 * Degree))), Tan(-42 * Degree));
}

TEST_F(QuantitiesTest, HyperbolicFunctions) {
  EXPECT_EQ(Sinh(0 * Radian), 0);
  EXPECT_EQ(Cosh(0 * Radian), 1);
  EXPECT_EQ(Tanh(0 * Radian), 0);
  // Limits:
  EXPECT_EQ(Sinh(20 * Radian), Cosh(20 * Radian));
  EXPECT_EQ(Tanh(20 * Radian), 1);
  EXPECT_EQ(Sinh(-20 * Radian), -Cosh(-20 * Radian));
  EXPECT_EQ(Tanh(-20 * Radian), -1);

  EXPECT_THAT(Sinh(2 * Radian) / Cosh(2 * Radian),
              AlmostEquals(Tanh(2 * Radian), 0, 1));
  EXPECT_THAT(ArcSinh(Sinh(-10 * Degree)),
              AlmostEquals(-10 * Degree, 0, 1));
  EXPECT_THAT(ArcCosh(Cosh(-10 * Degree)),
              AlmostEquals(10 * Degree, 19, 20));
  EXPECT_THAT(ArcTanh(Tanh(-10 * Degree)),
              AlmostEquals(-10 * Degree, 0, 1));
}

TEST_F(QuantitiesTest, ExpLogAndSqrt) {
  EXPECT_THAT(std::exp(std::log(2) / 2), AlmostEquals(Sqrt(2), 1));
  EXPECT_EQ(std::exp(std::log(Rood / Pow<2>(Foot)) / 2) * Foot, Sqrt(Rood));
}

TEST_F(QuantitiesDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Quantity message;
    message.set_dimensions(0x7C00);
    message.set_magnitude(1.0);
    Speed const speed_of_light = Speed::ReadFromMessage(message);
  }, "representation.*dimensions");
}

TEST_F(QuantitiesTest, SerializationSuccess) {
  serialization::Quantity message;
  SpeedOfLight.WriteToMessage(&message);
  EXPECT_EQ(0x7C01, message.dimensions());
  EXPECT_EQ(299792458.0, message.magnitude());
  Speed const speed_of_light = Speed::ReadFromMessage(message);
  EXPECT_EQ(SpeedOfLight, speed_of_light);
}

// This check verifies that setting a log handler causes the protobuf library to
// report its errors using glog.  It doesn't have much too do with quantities,
// except that it's a convenient protobuf for this test.
TEST_F(QuantitiesDeathTest, SerializationLogHandler) {
  EXPECT_DEATH({
    google::protobuf::SetLogHandler(
        [](google::protobuf::LogLevel const level,
           char const* const filename,
           int const line,
           std::string const& message) {
          LOG_AT_LEVEL(level) << "[" << filename << ":" << line << "] "
                              << message;
        });
    serialization::Quantity message;
    message.set_magnitude(1.0);
    message.CheckInitialized();
  }, "missing required fields");
}

}  // namespace quantities
}  // namespace principia
