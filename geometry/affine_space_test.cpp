#include "geometry/affine_space.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace geometry {

using geometry::Point;
using quantities::Temperature;
using si::Litre;
using si::Kelvin;
using testing::Eq;

class AffineSpaceTest : public testing::Test {
 protected:
};

Point<Temperature> const CelsiusZero(273.15 * Kelvin);
Point<Temperature> const AbsoluteZero(0 * Kelvin);

TEST_F(AffineSpaceTest, Comparisons) {
  EXPECT_TRUE(CelsiusZero == CelsiusZero);
  EXPECT_FALSE(CelsiusZero == AbsoluteZero);
  EXPECT_TRUE(CelsiusZero != AbsoluteZero);
  EXPECT_FALSE(CelsiusZero != CelsiusZero);
}

TEST_F(AffineSpaceTest, Operators) {
  Point<Temperature> const gallium_boiling_point =
      29.7646 * Kelvin + CelsiusZero;
  Point<Temperature> const water_boiling_point = CelsiusZero + 100 * Kelvin;
  Point<Temperature> const nitrogen_boiling_point =
      CelsiusZero - 195.795 * Kelvin;
  EXPECT_THAT(water_boiling_point - AbsoluteZero, Eq(373.15 * Kelvin));
  EXPECT_THAT(gallium_boiling_point - nitrogen_boiling_point,
              Eq(225.5596 * Kelvin));
}

TEST_F(AffineSpaceTest, Barycentres) {
  Point<Temperature> const T1 = 10 * Kelvin + CelsiusZero;
  Point<Temperature> const T2 = 40 * Kelvin + CelsiusZero;
  EXPECT_THAT(Barycentre(T1, 2 * Litre, T2, 1 * Litre),
              Eq(20 * Kelvin + CelsiusZero));
  EXPECT_THAT(Barycentre(T1, 1, T2, 1),
              Eq(25 * Kelvin + CelsiusZero));
}

}  // namespace geometry
}  // namespace principia
