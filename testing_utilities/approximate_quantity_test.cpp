
#include "testing_utilities/approximate_quantity.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_approximate_quantity {

TEST(ApproximateQuantityTest, Literals_⑴) {
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

TEST(ApproximateQuantityTest, Literals_⑵_⑼) {
  ApproximateQuantity<double> const l2 = 123.45_⑵;
  EXPECT_THAT(l2.min(), AlmostEquals(123.43, 0));
  EXPECT_THAT(l2.max(), AlmostEquals(123.47, 0));

  ApproximateQuantity<double> const l3 = 123.45_⑶;
  EXPECT_THAT(l3.min(), AlmostEquals(123.42, 0));
  EXPECT_THAT(l3.max(), AlmostEquals(123.48, 0));

  ApproximateQuantity<double> const l4 = 123.45_⑷;
  EXPECT_THAT(l4.min(), AlmostEquals(123.41, 0));
  EXPECT_THAT(l4.max(), AlmostEquals(123.49, 1));

  ApproximateQuantity<double> const l5 = 123.45_⑸;
  EXPECT_THAT(l5.min(), AlmostEquals(123.40, 0));
  EXPECT_THAT(l5.max(), AlmostEquals(123.50, 0));

  ApproximateQuantity<double> const l6 = 123.45_⑹;
  EXPECT_THAT(l6.min(), AlmostEquals(123.39, 0));
  EXPECT_THAT(l6.max(), AlmostEquals(123.51, 0));

  ApproximateQuantity<double> const l7 = 123.45_⑺;
  EXPECT_THAT(l7.min(), AlmostEquals(123.38, 1));
  EXPECT_THAT(l7.max(), AlmostEquals(123.52, 0));

  ApproximateQuantity<double> const l8 = 123.45_⑻;
  EXPECT_THAT(l8.min(), AlmostEquals(123.37, 0));
  EXPECT_THAT(l8.max(), AlmostEquals(123.53, 0));

  ApproximateQuantity<double> const l9 = 123.45_⑼;
  EXPECT_THAT(l9.min(), AlmostEquals(123.36, 0));
  EXPECT_THAT(l9.max(), AlmostEquals(123.54, 0));
}

}  // namespace internal_approximate_quantity
}  // namespace testing_utilities
}  // namespace principia
