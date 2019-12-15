
#include "geometry/rp2_point.hpp"

#include <limits>

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {
namespace internal_rp2_point {

using quantities::Infinity;
using quantities::Length;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;

class RP2PointTest : public ::testing::Test {
 protected:
  using Projective = Frame<serialization::Frame::TestTag,
                           Inertial,
                           Handedness::Right,
                           serialization::Frame::TEST>;
};

TEST_F(RP2PointTest, Basic) {
  RP2Point<Length, Projective> p1(1.0 * Metre, 2.0 * Metre, 3.0);
  RP2Point<Length, Projective> p2(2.0 * Metre, 4.0 * Metre, 6.0);
  RP2Point<Length, Projective> p3(2.0 * Metre, 4.0 * Metre, 5.0);
  RP2Point<Length, Projective> p4(1.0 * Metre, 2.0 * Metre, 0.0);
  RP2Point<Length, Projective> p5(2.0 * Metre, 4.0 * Metre, 0.0);
  RP2Point<Length, Projective> p6(0.0 * Metre, -4.0 * Metre, 0.0);
  RP2Point<Length, Projective> p7(0.0 * Metre, -0.0 * Metre, -0.0);

  // Basic equality.
  EXPECT_EQ(p1, p2);
  EXPECT_NE(p1, p3);
  EXPECT_EQ(p4, p5);
  EXPECT_NE(p4, p6);
  EXPECT_NE(p4, p7);

  // Infinities.
  EXPECT_FALSE(p1.is_at_infinity());
  EXPECT_FALSE(p2.is_at_infinity());
  EXPECT_FALSE(p3.is_at_infinity());
  EXPECT_TRUE(p4.is_at_infinity());
  EXPECT_TRUE(p5.is_at_infinity());
  EXPECT_TRUE(p6.is_at_infinity());
  EXPECT_TRUE(p7.is_at_infinity());

  // Euclidean coordinates.
  EXPECT_THAT(p1.x(), AlmostEquals(1.0 / 3.0 * Metre, 0));
  EXPECT_THAT(p1.y(), AlmostEquals(2.0 / 3.0 * Metre, 0));
  EXPECT_THAT(p2.x(), AlmostEquals(1.0 / 3.0 * Metre, 0));
  EXPECT_THAT(p2.y(), AlmostEquals(2.0 / 3.0 * Metre, 0));
  EXPECT_THAT(p3.x(), AlmostEquals(2.0 / 5.0 * Metre, 0));
  EXPECT_THAT(p3.y(), AlmostEquals(4.0 / 5.0 * Metre, 0));
  EXPECT_THAT(p4.x(), AlmostEquals(Infinity<Length>(), 0));
  EXPECT_THAT(p4.y(), AlmostEquals(Infinity<Length>(), 0));
  EXPECT_THAT(p5.x(), AlmostEquals(Infinity<Length>(), 0));
  EXPECT_THAT(p5.y(), AlmostEquals(Infinity<Length>(), 0));
  EXPECT_THAT(p6.x(), AlmostEquals(Infinity<Length>(), 0));
  EXPECT_THAT(p6.y(), AlmostEquals(-Infinity<Length>(), 0));
  EXPECT_THAT(p7.x(), AlmostEquals(-Infinity<Length>(), 0));
  EXPECT_THAT(p7.y(), AlmostEquals(Infinity<Length>(), 0));
}

}  // namespace internal_rp2_point
}  // namespace geometry
}  // namespace principia
