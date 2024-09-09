#include "numerics/angle_reduction.hpp"

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

template<typename T>
class AngleReductionTest : public testing::Test {};

TYPED_TEST_SUITE_P(AngleReductionTest);

// This test is not type-parameterized because the reduction algorithm only
// works for `Angle`.
TEST(AngleReductionTest, ReduceMinusπOver2ToπOver2) {
  Angle fractional_part;
  std::int64_t integer_part;

  ReduceAngle<-π / 2, π / 2>(Angle(1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(1 * Radian, 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-π / 2, π / 2>(Angle(-1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(-1 * Radian, 0));
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
      AlmostEquals(0.000030144353364053721297689416174085720 * Radian,
                   2221148));
  EXPECT_EQ(113, integer_part);
}

TYPED_TEST_P(AngleReductionTest, ReduceMinusπToπ) {
  using Angle = TypeParam;
  Angle fractional_part;
  std::int64_t integer_part;

  ReduceAngle<-π, π>(Angle(1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(1 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-π, π>(Angle(-1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(-1 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-π, π>(Angle(4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(4 * Radian) - two_π<Angle>, 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-π, π>(Angle(-4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(two_π<Angle> + Angle(-4 * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<-π, π>(Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(
      fractional_part,
      AlmostEquals(Angle(0x1.F9BD03091AD49p-15 * Radian) +
                       Angle(0x1.BA01B07B5D1EBp-71 * Radian),
                   4861461, 7230135));
  EXPECT_EQ(113, integer_part);
}

TYPED_TEST_P(AngleReductionTest, Reduce0To2π) {
  using Angle = TypeParam;
  Angle fractional_part;
  std::int64_t integer_part;

  ReduceAngle<0, 2 * π>(Angle(1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(1 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<0, 2 * π>(Angle(4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(4 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<0, 2 * π>(Angle(-1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(two_π<Angle> + Angle(-1 * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<0, 2 * π>(Angle(7 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(7 * Radian) - two_π<Angle>, 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-2 * π, 2 * π>(
      Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(
      fractional_part,
      AlmostEquals(Angle(0x1.F9BD03091AD49p-15 * Radian) +
                       Angle(0x1.BA01B07B5D1EBp-71 * Radian),
                   4861461, 7230135));
  EXPECT_EQ(113, integer_part);
}

TYPED_TEST_P(AngleReductionTest, ReduceMinus2πTo2π) {
  using Angle = TypeParam;
  Angle fractional_part;
  std::int64_t integer_part;

  ReduceAngle<-2 * π, 2 * π>(Angle(4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(4 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-2 * π, 2 * π>(Angle(-4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(-4 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<-2 * π, 2 * π>(Angle(7 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(7 * Radian) - two_π<Angle>, 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-2 * π, 2 * π>(Angle(-7 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(two_π<Angle> + Angle(-7 * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<-2 * π, 2 * π>(
      Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(
      fractional_part,
      AlmostEquals(Angle(0x1.F9BD03091AD49p-15 * Radian) +
                       Angle(0x1.BA01B07B5D1EBp-71 * Radian),
                   4861461, 7230135));
  EXPECT_EQ(113, integer_part);
}

REGISTER_TYPED_TEST_SUITE_P(AngleReductionTest,
                            ReduceMinusπToπ,
                            Reduce0To2π,
                            ReduceMinus2πTo2π);

using AngleTypes = ::testing::Types<Angle, DoublePrecision<Angle>>;
INSTANTIATE_TYPED_TEST_SUITE_P(AllAngleReductionTests,
                               AngleReductionTest,
                               AngleTypes);

}  // namespace internal
}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia
