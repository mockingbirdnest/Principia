#include "geometry/affine_space.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using geometry::Point;
using quantities::Temperature;
using si::Litre;
using si::Kelvin;
using testing::Eq;
using testing_utilities::AlmostEquals;

class AffineSpaceTest : public testing::Test {
 protected:
};

Point<Temperature> const CelsiusZero(273.15 * Kelvin);
Point<Temperature> const AbsoluteZero(0 * Kelvin);
Point<Temperature> const WaterBoilingPoint = CelsiusZero + 100 * Kelvin;
Point<Temperature> const GalliumMeltingPoint = 29.7646 * Kelvin + CelsiusZero;
Point<Temperature> const NitrogenBoilingPoint = CelsiusZero - 195.795 * Kelvin;

TEST_F(AffineSpaceTest, Comparisons) {
  EXPECT_TRUE(CelsiusZero == CelsiusZero);
  EXPECT_FALSE(CelsiusZero == AbsoluteZero);
  EXPECT_TRUE(CelsiusZero != AbsoluteZero);
  EXPECT_FALSE(CelsiusZero != CelsiusZero);
}

TEST_F(AffineSpaceTest, PlusMinus) {
  EXPECT_THAT(WaterBoilingPoint - AbsoluteZero, Eq(373.15 * Kelvin));
  EXPECT_THAT(GalliumMeltingPoint - NitrogenBoilingPoint,
              AlmostEquals(225.5596 * Kelvin));
}

TEST_F(AffineSpaceTest, AssignmentOperators) {
  Point<Temperature> accumulator = CelsiusZero;
  Point<Temperature> assignment_result;
  assignment_result = (accumulator += 100 * Kelvin);
  EXPECT_THAT(assignment_result, Eq(accumulator));
  EXPECT_THAT(accumulator, Eq(WaterBoilingPoint));
  assignment_result = (accumulator -= 100 * Kelvin);
  EXPECT_THAT(assignment_result, Eq(accumulator));
  EXPECT_THAT(accumulator, Eq(CelsiusZero));
  EXPECT_THAT((accumulator += 100 * Kelvin) -= 100 * Kelvin, Eq(CelsiusZero));
  EXPECT_THAT(accumulator, Eq(CelsiusZero));
}

TEST_F(AffineSpaceTest, Barycentres) {
  Point<Temperature> const t2 = 10 * Kelvin + CelsiusZero;
  Point<Temperature> const t1 = 40 * Kelvin + CelsiusZero;
  EXPECT_THAT(Barycentre(t2, 2 * Litre, t1, 1 * Litre),
              Eq(20 * Kelvin + CelsiusZero));
  EXPECT_THAT(Barycentre(t2, 1, t1, 1),
              Eq(25 * Kelvin + CelsiusZero));
}

}  // namespace geometry
}  // namespace principia
