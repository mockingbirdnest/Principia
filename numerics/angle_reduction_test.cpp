#include "numerics/angle_reduction.hpp"

#include <random>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

// The test is in the `internal` namespace to get visibility to `one_π` and
// `two_π`.
namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

template<typename Angle>
constexpr Angle one_π;

template<>
inline constexpr Angle one_π<Angle> = π * Radian;

template<>
inline constexpr DoublePrecision<Angle> one_π<DoublePrecision<Angle>> = []() {
  DoublePrecision<Angle> result;
  result.value = 0x1.921FB54442D18p1 * Radian;
  result.error = 0x1.1A62633145C07p-53 * Radian;
  return result;
}();

template<typename Angle>
constexpr Angle two_π;

template<>
inline constexpr Angle two_π<Angle> = 2 * π * Radian;

template<>
inline constexpr DoublePrecision<Angle> two_π<DoublePrecision<Angle>> = []() {
  DoublePrecision<Angle> result;
  result.value = 0x1.921FB54442D18p2 * Radian;
  result.error = 0x1.1A62633145C07p-52 * Radian;
  return result;
}();

class AngleReductionTest : public testing::Test {};

TEST_F(AngleReductionTest, PayneHanekMul97Examples) {
  {
    // [Mul97, Example 11, first angle].
    Angle const x = 0x1.8p200 * Radian;
    DoublePrecision<Angle> x_reduced;
    std::int64_t quadrant;
    PayneHanek<20>(x, x_reduced, quadrant);
    EXPECT_EQ(1, quadrant);
    // The last 20.4 bits of the result are incorrect.
    EXPECT_THAT(x_reduced,
                AlmostEquals(TwoSum(0x1.7F89C9C43D336p-1 * Radian,
                                    0x1.92CF93D957278p-56 * Radian),
                             1395298,
                             1395299));
  }
  {
    // [Mul97, Example 11, second angle].
    Angle const x = 6381956970095103.0 * 0x1.0p797 * Radian;
    DoublePrecision<Angle> x_reduced;
    std::int64_t quadrant;
    PayneHanek<61>(x, x_reduced, quadrant);
    EXPECT_EQ(1, quadrant);
    // The last 49.7 bits of the result are incorrect.
    EXPECT_THAT(x_reduced,
                AlmostEquals(TwoSum(0x1.14AE72E6BA22Fp-61 * Radian,
                                    -0x1.73EEF1477D90Ep-118 * Radian),
                             929348176455132,
                             929348176455133));
  }
}

TEST_F(AngleReductionTest, PayneHanekRandom) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> reduced_angle(-π / 4, π / 4);
  std::uniform_int_distribution<> quadrant(0, 3);
  std::uniform_int_distribution<> count(0, 10);
  std::uniform_int_distribution<> magnitude(0, 30);
  std::uniform_int_distribution<> sign(0, 1);
  for (std::int64_t i = 0; i < 1000; ++i) {
    // First, pick a reduced angle.
    DoublePrecision<Angle> const expected_reduced_angle(reduced_angle(random) *
                                                        Radian);
    std::int64_t const expected_quadrant = quadrant(random);
    // Then move it to a different quadrant by adding a multiple of  π / 2.
    DoublePrecision<Angle> x =
        expected_reduced_angle +
        expected_quadrant * Scale(0.5, one_π<DoublePrecision<Angle>>);
    // Now add a signed multiple of 2π.
    std::int64_t const loop_count = count(random);
    for (std::int64_t i = 0; i < loop_count; ++i) {
      auto const multiple_of_2π = Scale(std::scalbn(1.0, magnitude(random)),
                                        two_π<DoublePrecision<Angle>>);
      auto const signed_multiple_of_2π =
          (sign(random) * 2 - 1.0) * multiple_of_2π;
      x += signed_multiple_of_2π;
    }

    DoublePrecision<Angle> actual_reduced_angle;
    std::int64_t actual_quadrant;
    PayneHanek<61>(x.value, actual_reduced_angle, actual_quadrant);

    EXPECT_EQ(expected_quadrant, actual_quadrant);
    // We dropped `x.error` when calling `PayneHanek`, so we need to adjust our
    // expectations here.
    // The last 49.9 bits of the result may be incorrect, for
    // x = +1.68662971306440473e+09 rad|-5.42639714097227698e-08 rad, which has
    // a cancellation of 40.9 bits.  The fact that the error is so large is
    // because of the rounding errors on `x.error`, not because of the angle
    // reduction per se.
    EXPECT_THAT(actual_reduced_angle + x.error,
                AlmostEquals(expected_reduced_angle, 0, 1041371082701288))
        << "Expected reduced: " << expected_reduced_angle.value
        << " Reduction argument: " << x;
  }
}

// This test is not type-parameterized because the reduction algorithm only
// works for `Angle`.
TEST_F(AngleReductionTest, ReduceMinusπOver2ToπOver2) {
  Angle fractional_part;
  double integer_part;

  ReduceAngle<-π / 2, π / 2>(Angle(1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(1 * Radian, 1));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-π / 2, π / 2>(Angle(-1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(-1 * Radian, 1));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-π / 2, π / 2>(Angle(2 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals((2 - π) * Radian, 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-π / 2, π / 2>(Angle(-2 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals((π - 2) * Radian, 0));
  EXPECT_EQ(-1, integer_part);

  // 355/113 is hard, let's go shopping.
  ReduceAngle<-π / 2, π / 2>(
      Angle(355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(
      fractional_part,
      AlmostEquals(0.000030144353364053721297689416174085720 * Radian, 0));
  EXPECT_EQ(113, integer_part);
}

TEST_F(AngleReductionTest, ReduceMinusπToπ) {
  Angle fractional_part;
  double integer_part;

  ReduceAngle<-π, π>(Angle(1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(1 * Radian), 1));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-π, π>(Angle(-1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(-1 * Radian), 1));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-π, π>(Angle(3.5 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(3.5 * Radian) - two_π<Angle>, 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-π, π>(Angle(4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(4 * Radian) - two_π<Angle>, 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-π, π>(Angle(-4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(two_π<Angle> + Angle(-4 * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<-π, π>(Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(0x1.F9BD03091AD49p-15 * Radian) +
                               Angle(0x1.BA01B07B5D1EBp-71 * Radian),
                           0));
  EXPECT_EQ(113, integer_part);
}

TEST_F(AngleReductionTest, Reduce0To2π) {
  Angle fractional_part;
  double integer_part;

  ReduceAngle<0.0, 2 * π>(Angle(1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(1 * Radian), 1));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<0.0, 2 * π>(Angle(3.5 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(3.5 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<0.0, 2 * π>(Angle(4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(4 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<0.0, 2 * π>(Angle(-1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(two_π<Angle> + Angle(-1 * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<0.0, 2 * π>(Angle(-0.5 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(two_π<Angle> + Angle(-0.5 * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<0.0, 2 * π>(Angle(7 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(7 * Radian) - two_π<Angle>, 2));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<0.0, 2 * π>(
      Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(0x1.F9BD03091AD49p-15 * Radian) +
                               Angle(0x1.BA01B07B5D1EBp-71 * Radian),
                           0));
  EXPECT_EQ(113, integer_part);
}

}  // namespace internal
}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia
