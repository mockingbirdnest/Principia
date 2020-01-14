
#include "numerics/cbrt.hpp"

#include <cfenv>
#include <pmmintrin.h>

#include <limits>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/double_precision.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Truly;

#define EXPECT_SIGNALS(expression, exceptions)                            \
  do {                                                                    \
    std::feclearexcept(FE_ALL_EXCEPT);                                    \
    [[maybe_unused]] auto const volatile evaluated_signaling_expression = \
        (expression);                                                     \
    EXPECT_THAT(std::fetestexcept(FE_ALL_EXCEPT), Eq((exceptions)))       \
        << "while evaluating " #expression;                               \
  } while (false)

class CubeRootTest : public ::testing::Test {
 protected:
  CubeRootTest()
      : quiet_dead_beef_(FromBits(0x7FF8'0000'DEAD'BEEF)),
        signaling_dead_beef_(FromBits(0x7FF0'0000'DEAD'BEEF)) {}

  static std::uint64_t Bits(double const x) {
    return _mm_cvtsi128_si64(_mm_castpd_si128(_mm_set_sd(x)));
  }

  static double FromBits(std::uint64_t const x) {
    return _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128(x)));
  }

  static double ULP(double const x) {
    return FromBits(Bits(x) + 1) - x;
  }

  double const quiet_dead_beef_;
  // Microsoft's std::numeric_limits<double>::signaling_NaN() is actually quiet,
  // so we make our own.
  double const signaling_dead_beef_;
};

TEST_F(CubeRootTest, Cbrt2) {
  EXPECT_THAT(Cbrt(2), Eq(0x1.4'28A2'F98D'728Bp0));
  EXPECT_THAT(Cbrt(2) * Cbrt(2) * Cbrt(2), Eq(2));
}

TEST_F(CubeRootTest, Rescaling) {
  EXPECT_THAT(0x1p-340 * Cbrt(0x1p1021), Eq(Cbrt(2)));
  EXPECT_THAT(0x1p341 * Cbrt(0x1p-1022), Eq(Cbrt(2)));
  EXPECT_THAT(0x1p358 * Cbrt(0x1p-1073), Eq(Cbrt(2)));
}

TEST_F(CubeRootTest, SpecialValues) {
  EXPECT_THAT(Bits(Cbrt(0)), Eq(0));
  EXPECT_THAT(Bits(Cbrt(-0.0)), Eq(0x8000'0000'0000'0000));
  EXPECT_THAT(Cbrt(std::numeric_limits<double>::infinity()),
              Eq(std::numeric_limits<double>::infinity()));
  EXPECT_THAT(Cbrt(-std::numeric_limits<double>::infinity()),
              Eq(-std::numeric_limits<double>::infinity()));
  EXPECT_THAT(Cbrt(std::numeric_limits<double>::quiet_NaN()),
              Truly(&std::isnan<double>));
  EXPECT_THAT(Cbrt(-std::numeric_limits<double>::quiet_NaN()),
              Truly(&std::isnan<double>));
  // Preserve the payload of quiet NaNs as per IEEE 754-2008 6.2.
  EXPECT_THAT(Bits(Cbrt(quiet_dead_beef_)), Eq(Bits(quiet_dead_beef_)));
  EXPECT_THAT(Bits(Cbrt(signaling_dead_beef_)), Eq(Bits(quiet_dead_beef_)));
}

// This cube root is not correct as far as the inexact exception is concerned,
// but it should be OK otherwise (in other words, the only other exception it
// signals is the invalid operation exception on an sNaN operand).
TEST_F(CubeRootTest, Exceptions) {
  EXPECT_SIGNALS(Cbrt(2), FE_INEXACT);
  EXPECT_SIGNALS(Cbrt(8), FE_INEXACT);  // Erroneously.
  EXPECT_SIGNALS(Cbrt(signaling_dead_beef_), FE_INVALID);
  EXPECT_SIGNALS(Cbrt(quiet_dead_beef_), 0);
  EXPECT_SIGNALS(Cbrt(0), 0);
  EXPECT_SIGNALS(Cbrt(-0.0), 0);
  EXPECT_SIGNALS(Cbrt(std::numeric_limits<double>::infinity()), 0);
  EXPECT_SIGNALS(Cbrt(-std::numeric_limits<double>::infinity()), 0);
  EXPECT_SIGNALS(Cbrt(0x1p1021), FE_INEXACT);
  EXPECT_SIGNALS(Cbrt(0x1p-1022), FE_INEXACT);
  EXPECT_SIGNALS(Cbrt(0x1p-1073), FE_INEXACT);
}

TEST_F(CubeRootTest, BoundsOfTheRescalingRange) {
  EXPECT_THAT(Cbrt(0x1p-225), Eq(0x1p-75));
  EXPECT_THAT(Cbrt(0x1.0'0000'0000'0002p-225),
              Eq(0x1p-75 * Cbrt(0x1.0'0000'0000'0002p0)));
  EXPECT_THAT(Cbrt(0x1p237), Eq(0x1p79));
  EXPECT_THAT(Cbrt(0x1.F'FFFF'FFFF'FFFEp236),
              Eq(0x1p79 * Cbrt(0x1.F'FFFF'FFFF'FFFEp-1)));
}

TEST_F(CubeRootTest, Sign) {
  EXPECT_THAT(Cbrt(-2), Eq(-Cbrt(2)));
}

TEST_F(CubeRootTest, ParticularlyBadRounding) {
  constexpr double y = 0x1.14E35E87EA5DFp0;
  double const x = Cbrt(y);
  DoublePrecision<double> x² = TwoProduct(x, x);
  DoublePrecision<DoublePrecision<double>> x³ =
      TwoSum(TwoProduct(x².value, x), TwoProduct(x².error, x));
  CHECK_EQ(x³.error.error, 0) << x³;
  auto const long_y =
      DoublePrecision<DoublePrecision<double>>(DoublePrecision<double>(y));
  DoublePrecision<DoublePrecision<double>> numerator = long_y - x³;
  DoublePrecision<DoublePrecision<double>> denominator =
      TwoSum(TwoProduct(x².value, 3), TwoProduct(x².error, 3));
  double const x_error = numerator.value.value / denominator.value.value +
                         numerator.value.error / denominator.value.value;
  double const x_ulps = x_error / ULP(x);
  EXPECT_THAT(x_ulps, AllOf(Gt(0.5000551), Lt(0.5000552))) << x_ulps - 0.5;
}

}  // namespace numerics
}  // namespace principia
