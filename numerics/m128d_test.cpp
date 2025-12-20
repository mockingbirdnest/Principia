#include "numerics/m128d.hpp"

#include <cstdint>

#include "base/algebra.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/fma.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using namespace principia::base::_algebra;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_m128d;
using namespace principia::testing_utilities::_almost_equals;

class M128DTest : public testing::Test {};

TEST_F(M128DTest, Concepts) {
  // Not a field because 1 cannot be converted to M128D.
  static_assert(ring<M128D>);
}

TEST_F(M128DTest, Arithmetic) {
  M128D a(-5.0);
  M128D b(2.0);
  M128D c(3.0);

  EXPECT_THAT(static_cast<double>(-b), AlmostEquals(-2.0, 0));
  EXPECT_THAT(static_cast<double>(a + b), AlmostEquals(-3.0, 0));
  EXPECT_THAT(static_cast<double>(a - b), AlmostEquals(-7.0, 0));
  EXPECT_THAT(static_cast<double>(a * b), AlmostEquals(-10.0, 0));
  EXPECT_THAT(static_cast<double>(a / b), AlmostEquals(-2.5, 0));
  EXPECT_THAT(static_cast<double>(Abs(a)), AlmostEquals(5.0, 0));
  EXPECT_THAT(static_cast<double>(Abs(b)), AlmostEquals(2.0, 0));
  if constexpr (CanEmitFMAInstructions) {
    EXPECT_THAT(static_cast<double>(FusedMultiplyAdd(a, b, c)),
                AlmostEquals(-7.0, 0));
    EXPECT_THAT(static_cast<double>(FusedMultiplySubtract(a, b, c)),
                AlmostEquals(-13.0, 0));
    EXPECT_THAT(static_cast<double>(FusedNegatedMultiplyAdd(a, b, c)),
                AlmostEquals(13.0, 0));
    EXPECT_THAT(static_cast<double>(FusedNegatedMultiplySubtract(a, b, c)),
                AlmostEquals(7.0, 0));
  }

  a += b;
  EXPECT_THAT(static_cast<double>(a), AlmostEquals(-3.0, 0));
  a -= c;
  EXPECT_THAT(static_cast<double>(a), AlmostEquals(-6.0, 0));
  a *= b;
  EXPECT_THAT(static_cast<double>(a), AlmostEquals(-12.0, 0));
  a /= c;
  EXPECT_THAT(static_cast<double>(a), AlmostEquals(-4.0, 0));
}

TEST_F(M128DTest, Logical) {
  M128D a(5.0);
  M128D const sign_bit = M128D::MakeFromBits(0x8000'0000'0000'0000);
  EXPECT_EQ(0x8000'0000'0000'0000, sign_bit.Bits<std::uint64_t>());
  EXPECT_THAT(static_cast<double>(~a),
              AlmostEquals(-0.874999999999999888977697537484, 0));
  EXPECT_THAT(static_cast<double>(a & sign_bit), AlmostEquals(-0.0, 0));
  EXPECT_THAT(static_cast<double>(a | sign_bit), AlmostEquals(-5.0, 0));
  EXPECT_THAT(static_cast<double>(a ^ sign_bit), AlmostEquals(-5.0, 0));
}

TEST_F(M128DTest, Comparison) {
  M128D a(-5.0);
  M128D b(2.0);
  EXPECT_TRUE(a == a);
  EXPECT_TRUE(a == -5.0);
  EXPECT_TRUE(-5.0 == a);
  EXPECT_TRUE(a != b);
  EXPECT_TRUE(a != 2.0);
  EXPECT_TRUE(-5.0 != b);
  EXPECT_TRUE(a < b);
  EXPECT_TRUE(a < 2.0);
  EXPECT_TRUE(-5.0 < b);
  EXPECT_TRUE(a <= b);
  EXPECT_TRUE(a <= 2.0);
  EXPECT_TRUE(-5.0 <= b);
  EXPECT_TRUE(b >= a);
  EXPECT_TRUE(b >= -5.0);
  EXPECT_TRUE(2.0 >= a);
  EXPECT_TRUE(b > a);
  EXPECT_TRUE(b > -5.0);
  EXPECT_TRUE(2.0 > a);
}

}  // namespace numerics
}  // namespace principia
