
#include "numerics/cbrt.hpp"

#include <cfenv>
#include <pmmintrin.h>
#include <random>

#include <limits>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/double_precision.hpp"
#include "numerics/next.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {
namespace internal_cbrt {

using ::testing::AllOf;
using ::testing::AnyOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Ne;
using ::testing::Truly;

#define EXPECT_SIGNALS(expression, exceptions)                             \
  do {                                                                     \
    std::feclearexcept(FE_ALL_EXCEPT);                                     \
    [[maybe_unused]] auto const volatile evaluated_signaling_expression =  \
        (expression);                                                      \
    /* Ignore implementation-defined exceptions. */                        \
    EXPECT_THAT(std::fetestexcept(FE_DIVBYZERO | FE_INEXACT | FE_INVALID | \
                                  FE_OVERFLOW | FE_UNDERFLOW),             \
                Eq((exceptions)))                                          \
        << "while evaluating " #expression;                                \
  } while (false)

class CubeRootTest : public ::testing::Test {
 protected:
  CubeRootTest()
      : quiet_dead_beef_(FromBits(0x7FF8'0000'DEAD'BEEF)),
        signaling_dead_beef_(FromBits(0x7FF0'0000'DEAD'BEEF)) {}

  struct RoundedReal {
    double rounded_up;
    double rounded_down;
    double rounded_to_nearest;
  };

  // Uses a digit-by-digit algorithm to compute the correctly-rounded cube root
  // of y in all rounding modes.  y must lie in [1, 8[.
  static RoundedReal DigitByDigitCbrt(double const y) {
    CHECK_GE(y, 1);
    CHECK_LT(y, 8);
    double a = 1;
    double b = 0.5;
    for (int i = 1; i < 53; ++i, b /= 2) {
      a = CbrtOneBit(y, a, b) ? a + b : a;
    }
    bool const exact = TwoProduct(a, a).error == 0 &&
                       TwoProduct(a * a, a).error == 0 && a * a * a == y;
    RoundedReal result;
    result.rounded_down = a;
    result.rounded_to_nearest = CbrtOneBit(y, a, b) ? a + 2 * b : a;
    result.rounded_up = exact ? a : a + 2 * b;
    CHECK_LE(result.rounded_down, result.rounded_to_nearest);
    CHECK_LE(result.rounded_to_nearest, result.rounded_up);
    if (!exact) {
      CHECK_EQ(result.rounded_up, NextUp(result.rounded_down));
    }
    return result;
  }

  static void TestDigitByDigit(std::int64_t const values,
                               int const expected_3²ᴄZ5¹_misroundings,
                               int const expected_5²Z4¹FMA_misroundings) {
    // Properly testing the correction path is difficult because it is taken so
    // rarely; instead we test the bulk of it, namely CbrtOneBit, by using it to
    // compute the entire cube root and comparing against the actual
    // implementations.
    // This also allows us to test for faithfulness.

    std::mt19937_64 rng(1729);
    int method_3²ᴄZ5¹_misroundings = 0;
    int method_5²Z4¹FMA_misroundings = 0;
    for (std::int64_t i = 0; i < values; ++i) {
      double const y = FromBits(rng() % (Bits(8) - Bits(1)) + Bits(1));
      RoundedReal const cbrt_y = DigitByDigitCbrt(y);
      EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Faithful>(y),
                  AnyOf(cbrt_y.rounded_down, cbrt_y.rounded_up));
      EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(y),
                  Eq(cbrt_y.rounded_to_nearest));
      if (method_3²ᴄZ5¹::Cbrt<Rounding::Faithful>(y) !=
          cbrt_y.rounded_to_nearest) {
        ++method_3²ᴄZ5¹_misroundings;
      }
      if (CanEmitFMAInstructions) {
        EXPECT_THAT(method_5²Z4¹FMA::Cbrt<Rounding::Faithful>(y),
                    AnyOf(cbrt_y.rounded_down, cbrt_y.rounded_up));
        EXPECT_THAT(method_5²Z4¹FMA::Cbrt<Rounding::Correct>(y),
                    Eq(cbrt_y.rounded_to_nearest));
        if (method_5²Z4¹FMA::Cbrt<Rounding::Faithful>(y) !=
            cbrt_y.rounded_to_nearest) {
          ++method_5²Z4¹FMA_misroundings;
        }
      }
    }
    EXPECT_THAT(method_3²ᴄZ5¹_misroundings, Eq(expected_3²ᴄZ5¹_misroundings));
    if (CanEmitFMAInstructions) {
      EXPECT_THAT(method_5²Z4¹FMA_misroundings,
                  Eq(expected_5²Z4¹FMA_misroundings));
    }
  }

  static std::uint64_t Bits(double const x) {
    return _mm_cvtsi128_si64(_mm_castpd_si128(_mm_set_sd(x)));
  }

  static double FromBits(std::uint64_t const x) {
    return _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128(x)));
  }

  double const quiet_dead_beef_;
  // Microsoft's std::numeric_limits<double>::signaling_NaN() is actually quiet,
  // so we make our own.
  double const signaling_dead_beef_;
};

TEST_F(CubeRootTest, Rescaling) {
  EXPECT_THAT(0x1p-340 * method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(0x1p1021),
              Eq(Cbrt(2)));
  EXPECT_THAT(0x1p341 * method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(0x1p-1022),
              Eq(Cbrt(2)));
  EXPECT_THAT(0x1p358 * method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(0x1p-1073),
              Eq(Cbrt(2)));
  if (CanEmitFMAInstructions) {
    EXPECT_THAT(0x1p-340 * method_5²Z4¹FMA::Cbrt<Rounding::Correct>(0x1p1021),
                Eq(Cbrt(2)));
    EXPECT_THAT(0x1p341 * method_5²Z4¹FMA::Cbrt<Rounding::Correct>(0x1p-1022),
                Eq(Cbrt(2)));
    EXPECT_THAT(0x1p358 * method_5²Z4¹FMA::Cbrt<Rounding::Correct>(0x1p-1073),
                Eq(Cbrt(2)));
  }
}

TEST_F(CubeRootTest, DigitByDigit1e3) {
  TestDigitByDigit(/*values=*/1e3,
                   /*expected_3²ᴄZ5¹_misroundings=*/0,
                   /*expected_5²Z4¹FMA_misroundings=*/0);
}

#if !_DEBUG
TEST_F(CubeRootTest, DigitByDigit1e6) {
  TestDigitByDigit(/*values=*/1e6,
                   /*expected_3²ᴄZ5¹_misroundings=*/10,
                   /*expected_5²Z4¹FMA_misroundings=*/0);
}
#endif

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
  EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(0x1p-340),
              Eq(0x1p-114 * Cbrt(4)));
  EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(0x1.0'0000'0000'0002p-340),
              Eq(0x1p-114 * Cbrt(0x1.0'0000'0000'0002p2)));
  EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(0x1p341),
              Eq(0x1p113 * Cbrt(4)));
  EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(0x1.0'0000'0000'0002p341),
              Eq(0x1p113 * Cbrt(0x1.0'0000'0000'0002p2)));
  if (CanEmitFMAInstructions) {
    EXPECT_THAT(method_5²Z4¹FMA::Cbrt<Rounding::Correct>(0x1p-438),
                Eq(0x1p-146));
    EXPECT_THAT(
        method_5²Z4¹FMA::Cbrt<Rounding::Correct>(0x1.0'0000'0000'0002p-438),
        Eq(0x1p-146 * Cbrt(0x1.0'0000'0000'0002p0)));
    EXPECT_THAT(method_5²Z4¹FMA::Cbrt<Rounding::Correct>(0x1p438), Eq(0x1p146));
    EXPECT_THAT(
        method_5²Z4¹FMA::Cbrt<Rounding::Correct>(0x1.0'0000'0000'0002p438),
        Eq(0x1p146 * Cbrt(0x1.0'0000'0000'0002p0)));
  }
}

TEST_F(CubeRootTest, Sign) {
  EXPECT_THAT(Cbrt(-2), Eq(-Cbrt(2)));
}

TEST_F(CubeRootTest, ParticularlyDifficultRounding) {
  double const y = 0x1.1ACB56AEE37CEp1;
  // ∛y = 1.0100'1101'0110'1011'1110'0011'0110'1101'1010'1000'1001'1011'1001'
  //                       [1000'0000'0000'0000'0000'0000'0000'0000'0000'001…
  // with 37 0s after a 1 in the 54th bit.
  // Both faithful methods misround it.
  RoundedReal cbrt_y = DigitByDigitCbrt(y);
  EXPECT_THAT(cbrt_y.rounded_to_nearest,
              AllOf(Ne(cbrt_y.rounded_down), Eq(cbrt_y.rounded_up)));
  EXPECT_THAT(Cbrt(y), Eq(cbrt_y.rounded_to_nearest));
  EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Correct>(y),
              Eq(cbrt_y.rounded_to_nearest));
  EXPECT_THAT(method_3²ᴄZ5¹::Cbrt<Rounding::Faithful>(y),
              AllOf(Ne(cbrt_y.rounded_to_nearest), Eq(cbrt_y.rounded_down)));
  if (CanEmitFMAInstructions) {
    EXPECT_THAT(method_5²Z4¹FMA::Cbrt<Rounding::Correct>(y),
                Eq(cbrt_y.rounded_to_nearest));
    EXPECT_THAT(method_5²Z4¹FMA::Cbrt<Rounding::Faithful>(y),
                AllOf(Ne(cbrt_y.rounded_to_nearest), Eq(cbrt_y.rounded_down)));
  }
}

}  // namespace internal_cbrt
}  // namespace numerics
}  // namespace principia
