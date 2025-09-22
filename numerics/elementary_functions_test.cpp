#include "numerics/elementary_functions.hpp"

#include "base/cpuid.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "gtest/gtest.h"
#include "numerics/fma.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace numerics {

using ::testing::Eq;
using ::testing::Lt;
using namespace boost::multiprecision;
using namespace principia::base::_cpuid;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::quantities::_uk;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_vanishes_before;

class ElementaryFunctionsTest : public testing::Test {};

TEST_F(ElementaryFunctionsTest, FMA) {
  EXPECT_EQ(cpp_int(11), FusedMultiplyAdd(cpp_int(2), cpp_int(3), cpp_int(5)));
  EXPECT_EQ(cpp_rational(11, 2),
            FusedMultiplyAdd(
                cpp_rational(2, 1), cpp_rational(3, 2), cpp_rational(5, 2)));
  EXPECT_EQ(cpp_bin_float_50("11.0"),
            FusedMultiplyAdd(cpp_bin_float_50("2.0"),
                             cpp_bin_float_50("3.0"),
                             cpp_bin_float_50("5.0")));

  if (!CPUIDFeatureFlag::FMA.IsSet()) {
    GTEST_SKIP() << "Cannot test FMA on a machine without FMA";
  }
  EXPECT_EQ(11 * Coulomb,
            FusedMultiplyAdd(2 * Ampere, 3 * Second, 5 * Coulomb));
  EXPECT_EQ(11 * Radian, FusedMultiplyAdd(2.0, 3 * Radian, 5 * Radian));
  EXPECT_EQ(11.0, FusedMultiplyAdd(2.0, 3.0, 5.0));
}

TEST_F(ElementaryFunctionsTest, AbsoluteValue) {
  EXPECT_EQ(Abs(-1729), 1729);
  EXPECT_EQ(Abs(1729), 1729);
  EXPECT_EQ(Abs(-1729 * Metre), 1729 * Metre);
}

TEST_F(ElementaryFunctionsTest, DimensionlessExponentiation) {
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
  // This one calls `std::pow`.
  EXPECT_THAT(positivePower, AlmostEquals(Pow<4>(number), 0, 1));
  EXPECT_THAT(negativePower, AlmostEquals(Pow<-4>(number), 0, 1));
}

// The Greek letters cause a warning when stringified by the macros, because
// apparently Visual Studio doesn't encode strings in UTF-8 by default.
#pragma warning(disable: 4566)

TEST_F(ElementaryFunctionsTest, PhysicalConstants) {
  Length const lunar_distance = 384402 * Kilo(Metre);
  // By definition.
  EXPECT_THAT(1 / Pow<2>(SpeedOfLight),
              AlmostEquals(VacuumPermittivity * VacuumPermeability, 0));
  // The Keplerian approximation for the mass of the Sun
  // is fairly accurate.
  EXPECT_THAT(RelativeError(
                  4 * Pow<2>(π) * Pow<3>(AstronomicalUnit) / Pow<2>(JulianYear),
                  SolarGravitationalParameter),
              Lt(4e-5));
  EXPECT_THAT(RelativeError(1 * Parsec, 3.26156 * LightYear), Lt(2e-6));
  // The Keplerian approximation for the mass of the Earth
  // is pretty bad, but the error is still only 1%.
  EXPECT_THAT(RelativeError(4 * Pow<2>(π) * Pow<3>(lunar_distance) /
                                Pow<2>(27.321582 * Day),
                            TerrestrialGravitationalParameter),
              Lt(1e-2));
  EXPECT_THAT(RelativeError(1 * SolarGravitationalParameter,
                            1047 * JovianGravitationalParameter),
              Lt(6e-4));
  // Delambre & Méchain.
  EXPECT_THAT(RelativeError(TerrestrialGravitationalParameter /
                                Pow<2>(40 * Mega(Metre) / (2 * π)),
                            StandardGravity),
              Lt(4e-3));
  // Talleyrand.
  EXPECT_THAT(RelativeError(π * Sqrt(1 * Metre / StandardGravity), 1 * Second),
              Lt(4e-3));
}

#pragma warning(default: 4566)

TEST_F(ElementaryFunctionsTest, TrigonometricFunctions) {
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

TEST_F(ElementaryFunctionsTest, HyperbolicFunctions) {
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

TEST_F(ElementaryFunctionsTest, ExpLogAndRoots) {
  // The ULP distance is 1 if everything is correctly rounded.
  EXPECT_THAT(std::exp(std::log(2) / 2), AlmostEquals(Sqrt(2), 1));
  // The ULP distance is 0 if everything is correctly rounded.
  EXPECT_THAT(std::exp(std::log(2) / 3), AlmostEquals(Cbrt(2), 0));
  EXPECT_THAT(
      Sqrt(Rood),
      AlmostEquals(std::exp(std::log(Rood / Pow<2>(Foot)) / 2) * Foot, 0));
  EXPECT_THAT(
      Cbrt(Gallon),
      AlmostEquals(std::exp(std::log(Gallon / Pow<3>(Foot)) / 3) * Foot, 0, 1));
  auto const Lⵙ = SolarLuminosity;
  auto const Rⵙ = SolarRadius;
  auto const Tⵙ = Root<4>(Lⵙ / (4 * π * StefanBoltzmannConstant * Pow<2>(Rⵙ)));
  EXPECT_THAT(Tⵙ, IsNear(5772_(1) * Kelvin));
  EXPECT_THAT(Root(6, 100), AlmostEquals(Sqrt(Cbrt(100)), 0));
}

}  // namespace numerics
}  // namespace principia
