#include "geometry/plane.hpp"

#include "geometry/frame.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace geometry {

using ::testing::Eq;
using ::testing::UnorderedElementsAre;
using namespace principia::geometry::_plane;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

class PlaneTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag>;
};

TEST_F(PlaneTest, AmpèresLaw) {
  auto const I = Vector<Current, World>({2 * Ampere, -3 * Ampere, 5 * Ampere});
  auto const B =
      Bivector<MagneticFluxDensity, World>({-1 * Tesla, 0 * Tesla, 2 * Tesla});
  auto const plane_orthogonal_to_wire = Plane<World>::OrthogonalTo(I);
  auto const F = I * B;
  EXPECT_THAT(Projection(F, plane_orthogonal_to_wire), Eq(F));
}

TEST_F(PlaneTest, Rotation) {
  auto const Ω = AngularVelocity<World>(
      {2 * Radian / Second, -3 * Radian / Second, 5 * Radian / Second});
  auto const r = Displacement<World>({-1 * Metre, 0 * Metre, 2 * Metre});
  auto const plane_orthogonal_to_rotation_axis = Plane<World>::OrthogonalTo(Ω);
  auto const v = Ω * r;
  EXPECT_THAT(Projection(v, plane_orthogonal_to_rotation_axis), Eq(v));
}

TEST_F(PlaneTest, Normals) {
  auto const I = Vector<Current, World>({2 * Ampere, -3 * Ampere, 6 * Ampere});
  auto const plane_orthogonal_to_wire = Plane<World>::OrthogonalTo(I);
  EXPECT_THAT(
      plane_orthogonal_to_wire.UnitBinormals(),
      UnorderedElementsAre(
          Bivector<double, World>({2.0 / 7.0, -3.0 / 7.0, 6.0 / 7.0}),
          Bivector<double, World>({-2.0 / 7.0, 3.0 / 7.0, -6.0 / 7.0})));
  EXPECT_THAT(plane_orthogonal_to_wire.UnitNormals(),
              UnorderedElementsAre(
                  Vector<double, World>({2.0 / 7.0, -3.0 / 7.0, 6.0 / 7.0}),
                  Vector<double, World>({-2.0 / 7.0, 3.0 / 7.0, -6.0 / 7.0})));
}

}  // namespace geometry
}  // namespace principia
