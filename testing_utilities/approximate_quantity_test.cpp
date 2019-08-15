
#include "testing_utilities/approximate_quantity.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_approximate_quantity {

using quantities::Area;
using quantities::Frequency;
using quantities::Length;
using quantities::Speed;
using quantities::si::Metre;
using quantities::si::Second;

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

TEST(ApproximateQuantityTest, Units) {
  ApproximateQuantity<Length> const l1 = 123.45_⑴ * Metre;
  EXPECT_THAT(l1.min(), AlmostEquals(123.44 * Metre, 0));
  EXPECT_THAT(l1.max(), AlmostEquals(123.46 * Metre, 1));

  ApproximateQuantity<Area> const l2 = (123.45_⑴ * Metre) * Metre;
  EXPECT_THAT(l2.min(), AlmostEquals(123.44 * Metre * Metre, 0));
  EXPECT_THAT(l2.max(), AlmostEquals(123.46 * Metre * Metre, 1));

  ApproximateQuantity<Frequency> const l3 = 123.45_⑴ / Second;
  EXPECT_THAT(l3.min(), AlmostEquals(123.44 / Second, 0));
  EXPECT_THAT(l3.max(), AlmostEquals(123.46 / Second, 1));

  ApproximateQuantity<Speed> const l4 = 123.45_⑴ * Metre / Second;
  EXPECT_THAT(l4.min(), AlmostEquals(123.44 * Metre / Second, 0));
  EXPECT_THAT(l4.max(), AlmostEquals(123.46 * Metre / Second, 1));
}

TEST(ApproximateQuantityTest, DebugString) {
  EXPECT_EQ("[+1.23439999999999998e+02, +1.23460000000000008e+02]",
            (123.45_⑴).DebugString());
  EXPECT_EQ("[+1.23439999999999998e+02, +1.23460000000000008e+02] * "
            "+1.00000000000000000e+00 m",
            (123.45_⑴ * Metre).DebugString());
}

}  // namespace internal_approximate_quantity
}  // namespace testing_utilities
}  // namespace principia
