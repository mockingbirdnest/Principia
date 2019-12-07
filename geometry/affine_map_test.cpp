
#include "geometry/affine_map.hpp"

#include <limits>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/permutation.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace geometry {

using quantities::Length;
using quantities::si::Metre;
using quantities::si::Radian;
using testing::Contains;
using testing::Eq;
using testing::Lt;
using testing_utilities::AlmostEquals;
using testing_utilities::RelativeError;

class AffineMapTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  using Orth = OrthogonalMap<World, World>;
  using Perm = Permutation<World, World>;
  using Rot = Rotation<World, World>;
  using RigidTransformation = AffineMap<World, World, Length, Rotation>;

  AffineMapTest() {
    zero_     = Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre});
    forward_  = Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre});
    leftward_ = Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre});
    upward_   = Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre});

    origin_ = Position<World>();

    back_right_bottom_  = origin_ + Displacement<World>({3 * Metre,
                                                         4 * Metre,
                                                         5 * Metre});
    front_right_bottom_ = back_right_bottom_ + forward_;
    back_left_bottom_   = back_right_bottom_ + leftward_;
    front_left_bottom_  = back_left_bottom_ + forward_;
    back_right_top_     = back_right_bottom_ + upward_;
    front_right_top_    = back_right_top_ + forward_;
    back_left_top_      = back_right_top_ + leftward_;
    front_left_top_     = back_left_top_ + forward_;
    vertices_ = {back_left_bottom_, front_left_bottom_, back_right_bottom_,
                 front_right_bottom_, back_left_top_, front_left_top_,
                 back_right_top_, front_right_top_};
    originated_vertices_ = std::vector<Displacement<World>>(vertices_.size());
    for (std::size_t i = 0; i < vertices_.size(); ++i) {
      originated_vertices_[i] = vertices_[i] - origin_;
    }
  }

  Displacement<World> zero_;
  Displacement<World> forward_;
  Displacement<World> leftward_;
  Displacement<World> upward_;
  Position<World> origin_;
  Position<World> back_left_bottom_;
  Position<World> front_left_bottom_;
  Position<World> back_right_bottom_;
  Position<World> front_right_bottom_;
  Position<World> back_left_top_;
  Position<World> front_left_top_;
  Position<World> back_right_top_;
  Position<World> front_right_top_;
  std::vector<Position<World>> vertices_;
  std::vector<Displacement<World>> originated_vertices_;
};

TEST_F(AffineMapTest, Cube) {
  Rot const rotate_left(π / 2 * Radian,
                        Bivector<Length, World>(upward_.coordinates()));
  EXPECT_THAT(RelativeError(leftward_, rotate_left(forward_)),
              Lt(2 * std::numeric_limits<double>::epsilon()));
  RigidTransformation const map = RigidTransformation(back_right_bottom_,
                                                      front_right_bottom_,
                                                      rotate_left);
  EXPECT_THAT(map(back_right_bottom_) - origin_,
              Eq(front_right_bottom_ - origin_));
  EXPECT_THAT(map(front_left_top_) - origin_,
              AlmostEquals(back_left_top_ - origin_, 1));
  // Check that |map| is an isometry of the cube whose vertices are |vertices_|.
  for (auto const& point : vertices_) {
    EXPECT_THAT(originated_vertices_,
                Contains(AlmostEquals(map(point) - origin_, 0, 1)));
  }
  // Test that |map.Inverse() * map| acts as the identity on that cube.
  for (std::size_t i = 0; i < vertices_.size(); ++i) {
    EXPECT_THAT(originated_vertices_[i],
                AlmostEquals((map.Inverse() * map)
                                 (vertices_[i]) - origin_, 0));
  }
}

TEST_F(AffineMapTest, Serialization) {
  serialization::AffineMap message;
  Rot const rotate_left(π / 2 * Radian,
                        Bivector<Length, World>(upward_.coordinates()));
  RigidTransformation const map1 = RigidTransformation(back_right_bottom_,
                                                       front_right_bottom_,
                                                       rotate_left);
  map1.WriteToMessage(&message);
  EXPECT_TRUE(message.has_from_frame());
  EXPECT_TRUE(message.has_to_frame());
  EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
            message.to_frame().tag_type_fingerprint());
  EXPECT_EQ(message.from_frame().tag(),
            message.to_frame().tag());
  EXPECT_EQ(message.from_frame().is_inertial(),
            message.to_frame().is_inertial());
  EXPECT_TRUE(message.from_origin().has_multivector());
  EXPECT_TRUE(message.to_origin().has_multivector());
  EXPECT_TRUE(message.linear_map().HasExtension(
      serialization::Rotation::extension));
  RigidTransformation const map2 =
      RigidTransformation::ReadFromMessage(message);
  EXPECT_EQ(map1(back_right_bottom_) - origin_,
            map2(back_right_bottom_) - origin_);
  EXPECT_EQ(map1(front_left_top_) - origin_,
            map2(front_left_top_) - origin_);
}

TEST_F(AffineMapTest, Output) {
  Rot const rotate_left(π / 2 * Radian,
                        Bivector<Length, World>(upward_.coordinates()));
  RigidTransformation const map = RigidTransformation(back_right_bottom_,
                                                      front_right_bottom_,
                                                      rotate_left);
  std::cout << map << "\n";
}

}  // namespace geometry
}  // namespace principia
