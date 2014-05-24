#include "geometry/affine_map.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/permutation.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace geometry {

using quantities::Length;
using si::Metre;
using si::Radian;
using testing::Eq;
using testing_utilities::AlmostEquals;
using testing_utilities::RelativeError;

class AffineMapTest : public testing::Test {
 protected:
  struct World;
  typedef OrthogonalMap<World, World> Orth;
  typedef Permutation<World, World> Perm;
  typedef Rotation<World, World> Rot;
  typedef Vector<Length, World> Displacement;
  typedef Point<Displacement> Position;
  typedef AffineMap<World, World, Length, Rotation> RigidTransformation;

  void SetUp() override {
    zero_ = Displacement({0 * Metre, 0 * Metre, 0 * Metre});
    forward_ = Displacement({1 * Metre, 0 * Metre, 0 * Metre});
    rightward_ = Displacement({0 * Metre, 1 * Metre, 0 * Metre});
    upward_ = Displacement({0 * Metre, 0 * Metre, 1 * Metre});
    back_left_bottom_ = Position(Displacement({0 * Metre,
                                               0 * Metre,
                                               0 * Metre}));
    front_left_bottom_ = back_left_bottom_ + forward_;
    back_right_bottom_ = back_left_bottom_ + rightward_;
    front_right_bottom_ = back_right_bottom_ + forward_;
    back_left_top_ = back_left_bottom_ + upward_;
    front_left_top_ = back_left_top_ + forward_;
    back_right_top_ = back_left_top_ + rightward_;
    front_right_top_ = back_right_top_ + forward_;
  }

  Vector<Length, World> zero_;
  Vector<Length, World> forward_;
  Vector<Length, World> rightward_;
  Vector<Length, World> upward_;
  Position back_left_bottom_;
  Position front_left_bottom_;
  Position back_right_bottom_;
  Position front_right_bottom_;
  Position back_left_top_;
  Position front_left_top_;
  Position back_right_top_;
  Position front_right_top_;
};

TEST_F(AffineMapTest, Cube) {
  Rot rotate_right(-π / 2 * Radian, upward_);
  //EXPECT_THAT(RelativeError(rightward_, rotate_right(forward_)));
}

}  // namespace geometry
}  // namespace principia
