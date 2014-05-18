#include "geometry/affine_space.hpp"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace geometry {

using geometry::Point;
using quantities::Temperature;
using si::Kelvin;
using testing::Eq;

class AffineSpaceTest : public testing::Test {
 protected:
};

Point<Temperature> const CelsiusZero(273.15 * Kelvin);

TEST_F(AffineSpaceTest, WaterBoilingPoint) {
  Point<Temperature> const water_boiling_point = CelsiusZero + 100 * Kelvin;
  EXPECT_THAT(water_boiling_point - CelsiusZero, Eq(100 * Kelvin));
}

}  // namespace geometry
}  // namespace principia
