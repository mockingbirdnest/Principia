#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/permutation.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using quantities::Length;
using si::Metre;
using testing::Eq;
using testing_utilities::AlmostEquals;

class AffineMapTest : public testing::Test {
 protected:
  struct World;
  typedef OrthogonalMap<World, World> Orth;
  typedef Permutation<World, World> Perm;
  typedef Rotation<World, World> Rot;

  void SetUp() override {
    zero_ = Vector<Length, World>({0 * Metre, 0 * Metre, 0 * Metre});
    i_ = Vector<Length, World>({1 * Metre, 0 * Metre, 0 * Metre});
    j_ = Vector<Length, World>({0 * Metre, 1 * Metre, 0 * Metre});
    k_ = Vector<Length, World>({0 * Metre, 0 * Metre, 1 * Metre});
  }

  Vector<Length, World> zero_;
  Vector<Length, World> i_;
  Vector<Length, World> j_;
  Vector<Length, World> k_;
  Point<Length, World> bottom_left_back_;

};

TEST_F(AffineMapTest, Cube) {

}

}  // namespace geometry
}  // namespace principia
