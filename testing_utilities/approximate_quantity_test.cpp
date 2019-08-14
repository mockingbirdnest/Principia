
#include "testing_utilities/approximate_quantity.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_approximate_quantity {

TEST(ApproximateQuantityTest, Literals) {
  ApproximateQuantity<double> const l1 = 123.45_⑴;
  EXPECT_THAT(l1.min(), AlmostEquals(123.44, 0));
  EXPECT_THAT(l1.max(), AlmostEquals(123.46, 1));

  ApproximateQuantity<double> const l2 = 123.45e-3_⑴;
  EXPECT_THAT(l2.min(), AlmostEquals(123.44e-3, 1));
  EXPECT_THAT(l2.max(), AlmostEquals(123.46e-3, 0));

  ApproximateQuantity<double> const l3 = 123e3_⑴;
  EXPECT_THAT(l3.min(), AlmostEquals(122e3, 0));
  EXPECT_THAT(l3.max(), AlmostEquals(124e3, 0));

  ApproximateQuantity<double> const l4 = 0x1E3.45p0_⑴;
  EXPECT_THAT(l4.min(), AlmostEquals(0x1E3.44p0, 0));
  EXPECT_THAT(l4.max(), AlmostEquals(0x1E3.46p0, 0));

  ApproximateQuantity<double> const l5 = 123_⑴;
  EXPECT_THAT(l5.min(), AlmostEquals(122, 0));
  EXPECT_THAT(l5.max(), AlmostEquals(124, 0));
}

}  // namespace internal_approximate_quantity
}  // namespace testing_utilities
}  // namespace principia
