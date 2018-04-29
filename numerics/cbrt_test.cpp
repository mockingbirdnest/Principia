
#include "numerics/cbrt.hpp"

#include <pmmintrin.h>

#include <limits>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using ::testing::Eq;
using ::testing::Truly;

class CubeRootTest : public ::testing::Test {
 protected:
  static std::uint64_t bits(double const x) {
    return _mm_cvtsi128_si64(_mm_castpd_si128(_mm_set_sd(x)));
  }
};

TEST_F(CubeRootTest, Cbrt2) {
  EXPECT_THAT(cbrt(2), Eq(0x1.4'28A2'F98D'728Bp0));
  EXPECT_THAT(cbrt(2) * cbrt(2) * cbrt(2), Eq(2));
}

TEST_F(CubeRootTest, Rescaling) {
  EXPECT_THAT(0x1p-340 * cbrt(0x1p1021), Eq(cbrt(2)));
  EXPECT_THAT(0x1p341 * cbrt(0x1p-1022), Eq(cbrt(2)));
  EXPECT_THAT(0x1p358 * cbrt(0x1p-1073), Eq(cbrt(2)));
}

TEST_F(CubeRootTest, SpecialValues) {
  EXPECT_THAT(bits(cbrt(0)), Eq(0));
  EXPECT_THAT(bits(cbrt(-0.0)), Eq(0x8000'0000'0000'0000));
  EXPECT_THAT(cbrt(std::numeric_limits<double>::infinity()),
              Eq(std::numeric_limits<double>::infinity()));
  EXPECT_THAT(cbrt(-std::numeric_limits<double>::infinity()),
              Eq(-std::numeric_limits<double>::infinity()));
  EXPECT_THAT(cbrt(std::numeric_limits<double>::quiet_NaN()),
              Truly(&std::isnan<double>));
  EXPECT_THAT(cbrt(-std::numeric_limits<double>::quiet_NaN()),
              Truly(&std::isnan<double>));
}

TEST_F(CubeRootTest, BoundsOfTheRescalingRange) {
  EXPECT_THAT(cbrt(0x1p-225), Eq(0x1p-75));
  EXPECT_THAT(cbrt(0x1.0'0000'0000'0002p-225),
              Eq(0x1p-75 * cbrt(0x1.0'0000'0000'0002p0)));
  EXPECT_THAT(cbrt(0x1p237), Eq(0x1p79));
  EXPECT_THAT(cbrt(0x1.F'FFFF'FFFF'FFFFp236),
              Eq(0x1p79 * cbrt(0x1.F'FFFF'FFFF'FFFFp-1)));
}

TEST_F(CubeRootTest, Sign) {
  EXPECT_THAT(cbrt(-2), Eq(-cbrt(2)));
}

TEST_F(CubeRootTest, ParticularlyBadRounding) {
  EXPECT_THAT(cbrt(0x1.14E35E87EA5DFp0), Eq(0x1.06C80FCCA8E18p0));
}

}  // namespace numerics
}  // namespace principia
