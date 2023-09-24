#include "numerics/angle_reduction.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_angle_reduction;
using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

template<typename T>
class AngleReductionTest : public testing::Test {};

TYPED_TEST_SUITE_P(AngleReductionTest);

TEST(AngleReductionTest, SingleMinusπOver2ToπOver2) {
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
  ReduceAngle<-π / 2, π / 2>(Angle(355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals((355 - 113 * π) * Radian, 9451283));
  EXPECT_EQ(113, integer_part);
}

TYPED_TEST_P(AngleReductionTest, SingleMinusπToπ) {
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
  EXPECT_THAT(fractional_part, AlmostEquals(Angle((4 - 2 * π) * Radian), 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-π, π>(Angle(-4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle((2 * π - 4) * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<-π, π>(Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle((2 * 355 - 2 * 113 * π) * Radian), 0));
  EXPECT_EQ(113, integer_part);
}

TYPED_TEST_P(AngleReductionTest, Single0To2π) {
  using Angle = TypeParam;
  Angle fractional_part;
  std::int64_t integer_part;

  ReduceAngle<0, 2 * π>(Angle(4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(4 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<0, 2 * π>(Angle(4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle(4 * Radian), 0));
  EXPECT_EQ(0, integer_part);

  ReduceAngle<0, 2 * π>(Angle(-1 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle((2 * π - 1) * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<0, 2 * π>(Angle(7 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle((7 - 2 * π) * Radian), 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<0, 2 * π>(Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle((2 * 355 - 2 * 113 * π) * Radian), 0));
  EXPECT_EQ(113, integer_part);
}

TYPED_TEST_P(AngleReductionTest, SingleMinus2πTo2π) {
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
  EXPECT_THAT(fractional_part, AlmostEquals(Angle((7 - 2 * π) * Radian), 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-2 * π, 2 * π>(Angle(-7 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part, AlmostEquals(Angle((2 * π - 7) * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<-2 * π, 2 * π>(Angle(2 * 355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle((2 * 355 - 2 * 113 * π) * Radian), 0));
  EXPECT_EQ(113, integer_part);
}

REGISTER_TYPED_TEST_SUITE_P(AngleReductionTest,
                            SingleMinusπToπ,
                            Single0To2π,
                            SingleMinus2πTo2π);

using AngleTypes = ::testing::Types<Angle, DoublePrecision<Angle>>;
INSTANTIATE_TYPED_TEST_SUITE_P(AllAngleReductionTests,
                               AngleReductionTest,
                               AngleTypes);

}  // namespace numerics
}  // namespace principia
