// Start the file with a blank line or a comment to prevent linter confusion.
#include <string>

#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/astronomy.hpp"
#include "quantities/BIPM.hpp"
#include "quantities/constants.hpp"
#include "quantities/dimensionless.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/algebra.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/explicit_operators.hpp"
#include "testing_utilities/numerics.hpp"

using principia::astronomy::EarthMass;
using principia::astronomy::JulianYear;
using principia::astronomy::JupiterMass;
using principia::astronomy::LightYear;
using principia::astronomy::LunarDistance;
using principia::astronomy::Parsec;
using principia::astronomy::SolarMass;
using principia::constants::ElectronMass;
using principia::constants::GravitationalConstant;
using principia::constants::SpeedOfLight;
using principia::constants::StandardGravity;
using principia::constants::VacuumPermeability;
using principia::constants::VacuumPermittivity;
using principia::quantities::Abs;
using principia::quantities::Cos;
using principia::quantities::Dimensionless;
using principia::quantities::Exp;
using principia::quantities::Log;
using principia::quantities::Log10;
using principia::quantities::Log2;
using principia::quantities::Mass;
using principia::quantities::Product;
using principia::quantities::Sin;
using principia::quantities::Speed;
using principia::quantities::Sqrt;
using principia::si::Ampere;
using principia::si::AstronomicalUnit;
using principia::si::Candela;
using principia::si::Cycle;
using principia::si::Day;
using principia::si::Degree;
using principia::si::Hour;
using principia::si::Kelvin;
using principia::si::Kilogram;
using principia::si::Mega;
using principia::si::Metre;
using principia::si::Mole;
using principia::si::Radian;
using principia::si::Second;
using principia::si::Steradian;
using principia::testing_utilities::AbsoluteError;
using principia::testing_utilities::AlmostEquals;
using principia::testing_utilities::RelativeError;
using principia::testing_utilities::Times;
using principia::uk::Foot;
using principia::uk::Furlong;
using principia::uk::Mile;
using principia::uk::Rood;
using testing::Lt;

namespace principia {
namespace geometry {

class QuantitiesTest : public testing::Test {
 protected:
};

TEST_F(QuantitiesTest, AbsoluteValue) {
  EXPECT_EQ(Abs(-1729), 1729);
  EXPECT_EQ(Abs(1729), 1729);
}

TEST_F(QuantitiesTest, DimensionlessComparisons) {
  testing_utilities::TestOrder(Dimensionless(0), Dimensionless(1));
  testing_utilities::TestOrder(Dimensionless(-1), Dimensionless(0));
  testing_utilities::TestOrder(-e, e);
  testing_utilities::TestOrder(Dimensionless(3), π);
  testing_utilities::TestOrder(Dimensionless(42), Dimensionless(1729));
}

TEST_F(QuantitiesTest, DimensionfulComparisons) {
  testing_utilities::TestOrder(EarthMass, JupiterMass);
  testing_utilities::TestOrder(LightYear, Parsec);
  testing_utilities::TestOrder(-SpeedOfLight, SpeedOfLight);
  testing_utilities::TestOrder(SpeedOfLight * Day, LightYear);
}

TEST_F(QuantitiesTest, DimensionlessOperations) {
  Dimensionless const zero    = 0;
  Dimensionless const one     = 1;
  Dimensionless const taxi    = 1729;
  Dimensionless const answer  = 42;
  Dimensionless const heegner = 163;
  testing_utilities::TestField(
      zero, one, taxi, 4 * π / 3, heegner, answer, -e, 2);
}

TEST_F(QuantitiesTest, DimensionlfulOperations) {
  Dimensionless const zero = 0;
  Dimensionless const one = 1;
  Dimensionless const taxi = 1729;
  testing_utilities::TestVectorSpace(
      0 * Metre / Second, SpeedOfLight, 88 * Mile / Hour,
      -340.29 * Metre / Second, zero, one, -2 * π, taxi, 2);
  // Dimensionful multiplication is a tensor product, see [Tao 2012].
  testing_utilities::TestBilinearMap(
      Times<Product<Mass, Speed>, Mass, Speed>, SolarMass,
      ElectronMass, SpeedOfLight, 1 * Furlong / JulianYear, -e, 2);
}

TEST_F(QuantitiesTest, DimensionlessExponentiation) {
  Dimensionless const number   = π - 42;
  Dimensionless positivePowers = 1;
  Dimensionless negativePowers = 1;
  EXPECT_EQ(1, number.Pow<0>());
  for (int i = 1; i < 10; ++i) {
    positivePowers *= number;
    negativePowers /= number;
    EXPECT_THAT(number.Pow(i), AlmostEquals(positivePowers, i));
    EXPECT_THAT(number.Pow(-i), AlmostEquals(negativePowers, i));
  }
}

// The Greek letters cause a warning when stringified by the macros, because
// apparently Visual Studio doesn't encode strings in UTF-8 by default.
#pragma warning(disable: 4566)

TEST_F(QuantitiesTest, Formatting) {
  auto const allTheUnits = 1 * Metre * Kilogram * Second * Ampere * Kelvin /
                            (Mole * Candela * Cycle * Radian * Steradian);
  std::string const expected = std::string("1e+000 m kg s A K mol^-1") +
                                " cd^-1 cycle^-1 rad^-1 sr^-1";
  std::string const actual = ToString(allTheUnits, 0);
  EXPECT_EQ(expected, actual);
  std::string π16 = "3.1415926535897931e+000";
  EXPECT_EQ(ToString(π), π16);
}

TEST_F(QuantitiesTest, PhysicalConstants) {
  // By definition.
  EXPECT_THAT(1 / SpeedOfLight.Pow<2>(),
              AlmostEquals(VacuumPermittivity * VacuumPermeability, 2));
  // The Keplerian approximation for the mass of the Sun
  // is fairly accurate.
  EXPECT_THAT(RelativeError(
                  4 * π.Pow<2>() * AstronomicalUnit.Pow<3>() /
                      (GravitationalConstant * JulianYear.Pow<2>()),
                  SolarMass),
              Lt(4E-5));
  EXPECT_THAT(RelativeError(1 * Parsec, 3.26156 * LightYear), Lt(2E-6));
  // The Keplerian approximation for the mass of the Earth
  // is pretty bad, but the error is still only 1%.
  EXPECT_THAT(RelativeError(
                  4 * π.Pow<2>() * LunarDistance.Pow<3>() /
                      (GravitationalConstant * (27.321582 * Day).Pow<2>()),
                  EarthMass),
              Lt(1E-2));
  EXPECT_THAT(RelativeError(1 * SolarMass, 1047 * JupiterMass), Lt(4E-4));
  // Delambre & Méchain.
  EXPECT_THAT(RelativeError(
                  GravitationalConstant * EarthMass /
                      (40 * Mega(Metre) / (2 * π)).Pow<2>(),
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
  EXPECT_THAT(AbsoluteError(Cos(90 * Degree), 0), Lt(1E-16));
  EXPECT_EQ(Sin(90 * Degree), 1);
  EXPECT_EQ(Cos(180 * Degree), -1);
  EXPECT_THAT(AbsoluteError(Sin(180 * Degree), 0), Lt(1E-15));
  EXPECT_THAT(AbsoluteError(Cos(-90 * Degree), 0), Lt(1E-16));
  EXPECT_EQ(Sin(-90 * Degree), -1);
  for (int k = 1; k < 360; ++k) {
    // Don't test for multiples of 90 degrees as zeros lead to horrible
    // conditioning.
    if (k % 90 != 0) {
      EXPECT_THAT(Cos((90 - k) * Degree),
                  AlmostEquals(Sin(k * Degree), 50));
      EXPECT_THAT(Sin(k * Degree) / Cos(k * Degree),
                  AlmostEquals(Tan(k * Degree), 2));
      EXPECT_THAT(((k + 179) % 360 - 179) * Degree,
                  AlmostEquals(ArcTan(Sin(k * Degree), Cos(k * Degree)), 80));
      EXPECT_THAT(((k + 179) % 360 - 179) * Degree,
                  AlmostEquals(ArcTan(Sin(k * Degree) * AstronomicalUnit,
                                      Cos(k * Degree) * AstronomicalUnit), 80));
      EXPECT_THAT(Cos(ArcCos(Cos(k * Degree))),
                  AlmostEquals(Cos(k * Degree), 10));
      EXPECT_THAT(Sin(ArcSin(Sin(k * Degree))),
                  AlmostEquals(Sin(k * Degree), 1));
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

  EXPECT_EQ(Sinh(2 * Radian) / Cosh(2 * Radian), Tanh(2 * Radian));
  EXPECT_THAT(ArcSinh(Sinh(-10 * Degree)), AlmostEquals(-10 * Degree, 2));
  EXPECT_THAT(ArcCosh(Cosh(-10 * Degree)), AlmostEquals(10 * Degree, 20));
  EXPECT_THAT(ArcTanh(Tanh(-10 * Degree)), AlmostEquals(-10 * Degree, 1));
}

TEST_F(QuantitiesTest, ExpLogAndSqrt) {
  EXPECT_EQ(Exp(1), e);
  EXPECT_THAT(Exp(Log(4.2) + Log(1.729)), AlmostEquals(4.2 * 1.729, 1));
  EXPECT_EQ(Exp(Log(2) * Log2(1.729)), 1.729);
  EXPECT_EQ(Exp(Log(10) * Log10(1.729)), 1.729);
  EXPECT_THAT(Exp(Log(2) / 2), AlmostEquals(Sqrt(2), 1));
  EXPECT_EQ(Exp(Log(Rood / Foot.Pow<2>()) / 2) * Foot, Sqrt(Rood));
}

}  // namespace geometry
}  // namespace principia
