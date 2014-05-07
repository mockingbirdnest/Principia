#include "stdafx.hpp"

#include <CppUnitTest.h>
#include <stdio.h>

#include "Quantities/Astronomy.hpp"
#include "Quantities/BIPM.hpp"
#include "Quantities/Constants.hpp"
#include "Quantities/Dimensionless.hpp"
#include "Quantities/ElementaryFunctions.hpp"
#include "Quantities/Quantities.hpp"
#include "Quantities/SI.hpp"
#include "Quantities/UK.hpp"
#include "TestUtilities/Algebra.hpp"
#include "TestUtilities/GeometryComparisons.hpp"
#include "TestUtilities/QuantityComparisons.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace quantities {

using namespace astronomy;
using namespace constants;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace si;
using namespace test_utilities;
using namespace uk;

TEST_CLASS(QuantitiesTests) {
 public:
  TEST_METHOD(AbsoluteValue) {
    Assert::AreEqual(Abs(-1729), Dimensionless(1729));
    Assert::AreEqual(Abs(1729), Dimensionless(1729));
  }
  TEST_METHOD(DimensionlessComparisons) {
    TestOrder(Dimensionless(0), Dimensionless(1));
    TestOrder(Dimensionless(-1), Dimensionless(0));
    TestOrder(-e, e);
    TestOrder(Dimensionless(3), π);
    TestOrder(Dimensionless(42), Dimensionless(1729));
  }
  TEST_METHOD(DimensionfulComparisons) {
    TestOrder(EarthMass, JupiterMass);
    TestOrder(LightYear, Parsec);
    TestOrder(-SpeedOfLight, SpeedOfLight);
    TestOrder(SpeedOfLight * Day, LightYear);
  }
  TEST_METHOD(DimensionlessOperations) {
    Dimensionless const zero    = 0;
    Dimensionless const one     = 1;
    Dimensionless const taxi    = 1729;
    Dimensionless const answer  = 42;
    Dimensionless const heegner = 163;
    TestField(zero, one, taxi, 4 * π / 3, heegner, answer, -e);
  }
  TEST_METHOD(DimensionlfulOperations) {
    Dimensionless const zero = 0;
    Dimensionless const one = 1;
    Dimensionless const taxi = 1729;
    TestVectorSpace(0 * Metre / Second, SpeedOfLight, 88 * Mile / Hour,
                    -340.29 * Metre / Second, zero, one, -2 * π, taxi);
    auto multiply = [](Mass left, Speed right) { return left * right; };
    // Dimensionful multiplication is a tensor product, see [Tao 2012].
    TestBilinearMap(multiply, SolarMass, ElectronMass, SpeedOfLight,
                    1 * Furlong / JulianYear, -e);
  }
  TEST_METHOD(DimensionlessExponentiation) {
    Dimensionless const number   = π - 42;
    Dimensionless positivePowers = 1;
    Dimensionless negativePowers = 1;
    AssertEqual(Dimensionless(1), number.Pow<0>());
    for (int i = 1; i < 10; ++i) {
      positivePowers *= number;
      negativePowers /= number;
      AssertEqual(number.Pow(i), positivePowers);
      AssertEqual(number.Pow(-i), negativePowers);
    }
  }
  TEST_METHOD(Formatting) {
    auto const allTheUnits = 1 * Metre * Kilogram * Second * Ampere * Kelvin /
                             (Mole * Candela * Cycle * Radian * Steradian);
    std::wstring const expected = std::wstring(L"1e+000 m^1 kg^1 s^1") + 
                                  L" A^1 K^1 mol^-1 cd^-1 cycle^-1 rad^-1" +
                                  L" sr^-1";
    std::wstring const actual = ToString(allTheUnits, 0);
    AssertTrue(actual.compare(expected) == 0);
    std::wstring π16 = L"3.1415926535897931e+000";
    AssertTrue(ToString(π).compare(π16) == 0);
  }
  TEST_METHOD(PhysicalConstants) {
    AssertEqual(1 / SpeedOfLight.Pow<2>(),
                VacuumPermittivity * VacuumPermeability);
    // The Keplerian approximation for the mass of the Sun
    // is fairly accurate.
    AssertEqual(4 * π.Pow<2>() * AstronomicalUnit.Pow<3>() /
                    (GravitationalConstant * JulianYear.Pow<2>()),
                SolarMass, 1e-4);
    AssertEqual(1 * Parsec, 3.26156 * LightYear, 1e-5);
    // The Keplerian approximation for the mass of the Earth
    // is pretty bad, but the error is still only 1%.
    AssertEqual(4 * π.Pow<2>() * LunarDistance.Pow<3>() /
                    (GravitationalConstant * (27.321582 * Day).Pow<2>()),
                EarthMass, 1e-2);
    AssertEqual(1 * SolarMass, 1047 * JupiterMass, 1e-3);
    // Delambre & Méchain.
    AssertEqual(GravitationalConstant * EarthMass /
                    (40 * Mega(Metre) / (2 * π)).Pow<2>(),
                StandardGravity, 1e-2);
    // Talleyrand.
    AssertEqual(π * Sqrt(1 * Metre / StandardGravity),
                1 * Second, 1e-2);
  }
  TEST_METHOD(TrigonometricFunctions) {
    AssertEqual(Cos(0 * Degree), 1);
    AssertEqual(Sin(0 * Degree), 0);
    AssertEqual(Cos(90 * Degree), 0);
    AssertEqual(Sin(90 * Degree), 1);
    AssertEqual(Cos(180 * Degree), -1);
    AssertEqual(Sin(180 * Degree), 0);
    AssertEqual(Cos(-90 * Degree), 0);
    AssertEqual(Sin(-90 * Degree), -1);
    for (int k = 0; k < 360; ++k) {
      AssertEqualAbsolute(Cos((90 - k) * Degree), Sin(k * Degree));
      AssertEqual(Sin(k * Degree) / Cos(k * Degree), Tan(k * Degree));
      AssertEqual(((k + 179) % 360 - 179) * Degree,
                  ArcTan(Sin(k * Degree), Cos(k * Degree)),
                  1e-13);
      AssertEqual(((k + 179) % 360 - 179) * Degree,
                  ArcTan(Sin(k * Degree) * AstronomicalUnit,
                  Cos(k * Degree) * AstronomicalUnit),
                  1e-13);
      AssertEqualAbsolute(Cos(ArcCos(Cos(k * Degree))), Cos(k * Degree));
      AssertEqualAbsolute(Sin(ArcSin(Sin(k * Degree))), Sin(k * Degree));
      AssertEqual(Sin(k * Degree).Pow(2), (1 - Cos(2 * k * Degree)) / 2, 1e-13);
      AssertEqual(Cos(k * Degree).Pow(2), (1 + Cos(2 * k * Degree)) / 2, 1e-13);
    }
    // Horribly conditioned near 0, so not in the loop above.
    AssertEqual(Tan(ArcTan(Tan(-42 * Degree))), Tan(-42 * Degree));
  }
  TEST_METHOD(HyperbolicFunctions) {
    AssertEqual(Sinh(0 * Radian), 0);
    AssertEqual(Cosh(0 * Radian), 1);
    AssertEqual(Tanh(0 * Radian), 0);
    // Limits:
    AssertEqual(Sinh(20 * Radian), Cosh(20 * Radian));
    AssertEqual(Tanh(20 * Radian), 1);
    AssertEqual(Sinh(-20 * Radian), -Cosh(-20 * Radian));
    AssertEqual(Tanh(-20 * Radian), -1);
    
    AssertEqual(Sinh(2 * Radian) / Cosh(2 * Radian), Tanh(2 * Radian));
    AssertEqual(ArcSinh(Sinh(-10 * Degree)), -10 * Degree);
    AssertEqual(ArcCosh(Cosh(-10 * Degree)), 10 * Degree, 1e-14);
    AssertEqual(ArcTanh(Tanh(-10 * Degree)), -10 * Degree);
  }
  TEST_METHOD(ExpLogAndSqrt) {
    AssertEqual(Exp(1), e);
    AssertEqual(Exp(Log(4.2) + Log(1.729)), 4.2 * 1.729, 1e-14);
    AssertEqual(Exp(Log(2) * Log2(1.729)), 1.729);
    AssertEqual(Exp(Log(10) * Log10(1.729)), 1.729);
    AssertEqual(Exp(Log(2) / 2), Sqrt(2));
    AssertEqual(Exp(Log(Rood / Foot.Pow<2>()) / 2) * Foot, Sqrt(Rood));
  }
};

}  // namespace quantities
}  // namespace principia
