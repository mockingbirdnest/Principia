
#include "geometry/r3_projective.hpp"

#include <limits>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {
namespace internal_r3_projective {

using quantities::Infinity;
using quantities::Length;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;

class R3ProjectiveTest : public ::testing::Test {};

TEST_F(R3ProjectiveTest, Basic) {
  R3Projective<Length> p1({1.0 * Metre, 2.0 * Metre, 3.0});
  R3Projective<Length> p2({2.0 * Metre, 4.0 * Metre, 6.0});
  R3Projective<Length> p3({2.0 * Metre, 4.0 * Metre, 5.0});
  R3Projective<Length> p4({1.0 * Metre, 2.0 * Metre, 0.0});
  R3Projective<Length> p5({2.0 * Metre, 4.0 * Metre, 0.0});
  R3Projective<Length> p6({0.0 * Metre, -4.0 * Metre, 0.0});
  auto b = p1 == p2;

  // Basic equality.
  EXPECT_EQ(p1, p2);
  EXPECT_NE(p1, p3);
  EXPECT_EQ(p4, p5);
  EXPECT_NE(p4, p6);

  // Infinities.
  EXPECT_FALSE(p1.is_at_infinity());
  EXPECT_FALSE(p2.is_at_infinity());
  EXPECT_FALSE(p3.is_at_infinity());
  EXPECT_TRUE(p4.is_at_infinity());
  EXPECT_TRUE(p5.is_at_infinity());
  EXPECT_TRUE(p6.is_at_infinity());
  EXPECT_THAT(p4.point_at_infinity(), AlmostEquals(2.0, 0));
  EXPECT_THAT(p5.point_at_infinity(), AlmostEquals(2.0, 0));
  EXPECT_THAT(p6.point_at_infinity(),
              AlmostEquals(std::numeric_limits<double>::infinity(), 0));

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
}

}  // namespace internal_r3_projective
}  // namespace geometry
}  // namespace principia
