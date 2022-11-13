#include "numerics/next.hpp"

#include <limits>

#include "geometry/sign.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace numerics {

using geometry::Sign;
using ::testing::Eq;

class NextTest : public ::testing::Test {};

TEST_F(NextTest, ExtendedRealLine) {
  static_assert(NextUp(-std::numeric_limits<double>::infinity()) ==
                -0x1.FFFFFFFFFFFFFp+1023);
  static_assert(NextUp(-0x1.FFFFFFFFFFFFFp+1023) == -0x1.FFFFFFFFFFFFEp+1023);
  // [...]
  static_assert(NextUp(-0x1.0000000000001p+1023) == -0x1.0000000000000p+1023);
  // Exponent change between negative normal numbers below -1.
  static_assert(NextUp(-0x1.0000000000000p+1023) == -0x1.FFFFFFFFFFFFFp+1022);
  static_assert(NextUp(-0x1.FFFFFFFFFFFFFp+1022) == -0x1.FFFFFFFFFFFFEp+1022);
  // [...]
  static_assert(NextUp(-0x1.0000000000001p-1021) == -0x1.0000000000000p-1021);
  // Exponent change between negative normal numbers above -1.
  static_assert(NextUp(-0x1.0000000000000p-1021) == -0x1.FFFFFFFFFFFFFp-1022);
  static_assert(NextUp(-0x1.FFFFFFFFFFFFFp-1022) == -0x1.FFFFFFFFFFFFEp-1022);
  // [...]
  static_assert(NextUp(-0x1.0000000000002p-1022) == -0x1.0000000000001p-1022);
  static_assert(NextUp(-0x1.0000000000001p-1022) == -0x1.0000000000000p-1022);
  // Normal-to-subnormal transition.
  static_assert(NextUp(-0x1.0000000000000p-1022) == -0x0.FFFFFFFFFFFFFp-1022);
  static_assert(NextUp(-0x0.FFFFFFFFFFFFFp-1022) == -0x0.FFFFFFFFFFFFEp-1022);
  // [...]
  static_assert(NextUp(-0x0.0000000000002p-1022) == -0x0.0000000000001p-1022);
  static_assert(NextUp(-0x0.0000000000001p-1022) == 0);
  // That was a -0, but we cannot check that in constexpr code until we have
  // std::bit_cast in C++20.
  EXPECT_THAT(Sign(NextUp(-0x0.0000000000001p-1022)), Eq(Sign::Negative()));
  // NextUp(0), being the least number greater than 0, does not depend on the
  // sign of 0 (contrast with std::nextafter).
  static_assert(NextUp(-0.0) == 0x0.0000000000001p-1022);
  static_assert(NextUp(+0.0) == 0x0.0000000000001p-1022);
  static_assert(NextUp(0x0.0000000000001p-1022) == 0x0.0000000000002p-1022);
  // [...]
  static_assert(NextUp(0x0.FFFFFFFFFFFFEp-1022) == 0x0.FFFFFFFFFFFFFp-1022);
  // Subnormal-to-normal transition.
  static_assert(NextUp(0x0.FFFFFFFFFFFFFp-1022) == 0x1.0000000000000p-1022);
  static_assert(NextUp(0x1.0000000000000p-1022) == 0x1.0000000000001p-1022);
  static_assert(NextUp(0x1.0000000000001p-1022) == 0x1.0000000000002p-1022);
  // [...]
  static_assert(NextUp(0x1.FFFFFFFFFFFFEp-1022) == 0x1.FFFFFFFFFFFFFp-1022);
  // Exponent change between p positive numbers below 1.
  static_assert(NextUp(0x1.FFFFFFFFFFFFFp-1022) == 0x1.0000000000000p-1021);
  static_assert(NextUp(0x1.0000000000000p-1021) == 0x1.0000000000001p-1021);
  // [...]
  static_assert(NextUp(0x1.FFFFFFFFFFFFEp+1022) == 0x1.FFFFFFFFFFFFFp+1022);
  // Exponent change between positive normal numbers above 1.
  static_assert(NextUp(0x1.FFFFFFFFFFFFFp+1022) == 0x1.0000000000000p+1023);
  static_assert(NextUp(0x1.0000000000000p+1023) == 0x1.0000000000001p+1023);
  // [...]
  static_assert(NextUp(0x1.FFFFFFFFFFFFEp+1023) == 0x1.FFFFFFFFFFFFFp+1023);
  static_assert(NextUp(0x1.FFFFFFFFFFFFFp+1023) ==
                std::numeric_limits<double>::infinity());
  // The train stops here.
  static_assert(NextUp(std::numeric_limits<double>::infinity()) ==
                std::numeric_limits<double>::infinity());
}

TEST_F(NextTest, MuchAdoAboutNothingsSignBit) {
  // Constructing negative 0s by going below 0 and back up to it.
  static_assert(NextUp(NextDown(0.0)) == 0);
  EXPECT_THAT(Sign(NextUp(NextDown(0.0))), Eq(Sign::Negative()));
}

}  // namespace numerics
}  // namespace principia
