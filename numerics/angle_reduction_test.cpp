#include "numerics/angle_reduction.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

// The test is in the |internal| namespace to get visibility to |one_π| and
// |two_π|.
namespace principia {
namespace numerics {
namespace _angle_reduction {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

template<typename T>
class AngleReductionTest : public testing::Test {
 protected:
  static inline T const four_π = two_π<T> + two_π<T>;
  static inline T const eight_π = four_π + four_π;
  static inline T const sixteen_π = eight_π + eight_π;
  static inline T const thirtytwo_π = sixteen_π + sixteen_π;
  static inline T const sixtyfour_π = thirtytwo_π + thirtytwo_π;
  static inline T const hundredandtwelve_π =
      sixtyfour_π + thirtytwo_π + sixteen_π;
  static inline T const hundredandfourteen_π =
      sixtyfour_π + thirtytwo_π + sixteen_π + two_π<T>;
  static inline T const twohundredandtwentysix_π =
      hundredandtwelve_π + hundredandfourteen_π;
};

TYPED_TEST_SUITE_P(AngleReductionTest);

// This test is not type-parameterized because the reduction algorithm only
// works for |Angle|.
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
  ReduceAngle<-π / 2, π / 2>(
      Angle(355 * Radian), fractional_part, integer_part);
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
  EXPECT_THAT(fractional_part,
              AlmostEquals(Angle(4 * Radian) - two_π<Angle>, 0));
  EXPECT_EQ(1, integer_part);

  ReduceAngle<-π, π>(Angle(-4 * Radian), fractional_part, integer_part);
  EXPECT_THAT(fractional_part,
              AlmostEquals(two_π<Angle> + Angle(-4 * Radian), 0));
  EXPECT_EQ(-1, integer_part);

  ReduceAngle<-π, π>(Angle(355 * Radian), fractional_part, integer_part);
  EXPECT_THAT(
      fractional_part,
      AlmostEquals(
          Angle(355 * Radian) - TestFixture::hundredandfourteen_π, 14, 48));
  EXPECT_EQ(57, integer_part);
}

TYPED_TEST_P(AngleReductionTest, Single0To2π) {
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
      AlmostEquals(
          Angle(2 * 355 * Radian) - TestFixture::twohundredandtwentysix_π, 0));
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
      AlmostEquals(
          Angle(2 * 355 * Radian) - TestFixture::twohundredandtwentysix_π, 0));
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

}  // namespace internal
}  // namespace _angle_reduction
}  // namespace numerics
}  // namespace principia
