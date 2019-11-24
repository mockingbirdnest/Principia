
#include "geometry/orthogonal_map.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/rotation.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {
namespace internal_orthogonal_map {

using quantities::si::Degree;
using quantities::si::Metre;
using testing::Eq;
using testing_utilities::AlmostEquals;

class OrthogonalMapTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;
  using Orth = OrthogonalMap<World, World>;
  using Rot = Rotation<World, World>;

  OrthogonalMapTest()
      : vector_(Vector<quantities::Length, World>(
            R3Element<quantities::Length>(
                1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        bivector_(Bivector<quantities::Length, World>(
            R3Element<quantities::Length>(
                1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        trivector_(Trivector<quantities::Length, World>(4.0 * Metre)),
        orthogonal_a_(Orth(Sign::Negative(),
                           Rot(120 * Degree,
                               Bivector<double, World>({1, 1, 1})))),
        orthogonal_b_(Orth(Sign::Positive(),
                           Rot(90 * Degree,
                               Bivector<double, World>({1, 0, 0})))),
        orthogonal_c_(Orth(Sign::Negative(),
                           Rot(90 * Degree,
                               Bivector<double, World>({1, 0, 0})))) {}

  Vector<quantities::Length, World> vector_;
  Bivector<quantities::Length, World> bivector_;
  Trivector<quantities::Length, World> trivector_;
  Orth orthogonal_a_;
  Orth orthogonal_b_;
  Orth orthogonal_c_;
};

using OrthogonalMapDeathTest = OrthogonalMapTest;

TEST_F(OrthogonalMapTest, Identity) {
  EXPECT_THAT(vector_, Eq(Orth::Identity()(vector_)));
  EXPECT_THAT(bivector_, Eq(Orth::Identity()(bivector_)));
  EXPECT_THAT(trivector_, Eq(Orth::Identity()(trivector_)));
}

TEST_F(OrthogonalMapTest, AppliedToVector) {
  EXPECT_THAT(orthogonal_a_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(-3.0 * Metre,
                                                -1.0 * Metre,
                                                -2.0 * Metre)), 4));
  EXPECT_THAT(orthogonal_b_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1, 2));
}

TEST_F(OrthogonalMapTest, AppliedToBivector) {
  EXPECT_THAT(orthogonal_a_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
  EXPECT_THAT(orthogonal_b_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1, 2));
}

TEST_F(OrthogonalMapTest, AppliedToTrivector) {
  EXPECT_THAT(orthogonal_a_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  -4.0 * Metre), 0));
  EXPECT_THAT(orthogonal_b_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
}

TEST_F(OrthogonalMapTest, Determinant) {
  EXPECT_TRUE(orthogonal_a_.Determinant().is_negative());
  EXPECT_TRUE(orthogonal_b_.Determinant().is_positive());
  EXPECT_TRUE(orthogonal_c_.Determinant().is_negative());
}

TEST_F(OrthogonalMapTest, Inverse) {
  EXPECT_THAT(orthogonal_a_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(-2.0 * Metre,
                                                -3.0 * Metre,
                                                -1.0 * Metre)), 2));
  EXPECT_THAT(orthogonal_b_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre)), 1, 2));
}

TEST_F(OrthogonalMapTest, Composition) {
  Orth const orthogonal_ac = orthogonal_a_ * orthogonal_c_;
  EXPECT_THAT(orthogonal_ac(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                1.0 * Metre,
                                                -3.0 * Metre)), 4, 6));
  EXPECT_TRUE((orthogonal_a_ * orthogonal_b_).Determinant().is_negative());
  EXPECT_TRUE((orthogonal_a_ * orthogonal_c_).Determinant().is_positive());
  EXPECT_TRUE((orthogonal_b_ * orthogonal_c_).Determinant().is_negative());
}

TEST_F(OrthogonalMapDeathTest, SerializationError) {
  Identity<World, World> id;
  EXPECT_DEATH({
    serialization::LinearMap message;
    id.WriteToMessage(&message);
    Orth const o = Orth::ReadFromMessage(message);
  }, "HasExtension.*OrthogonalMap");
}

TEST_F(OrthogonalMapTest, SerializationSuccess) {
  serialization::LinearMap message;
  orthogonal_a_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_from_frame());
  EXPECT_TRUE(message.has_to_frame());
  EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
            message.to_frame().tag_type_fingerprint());
  EXPECT_EQ(message.from_frame().tag(),
            message.to_frame().tag());
  EXPECT_EQ(message.from_frame().is_inertial(),
            message.to_frame().is_inertial());
  EXPECT_TRUE(message.HasExtension(
      serialization::OrthogonalMap::extension));
  serialization::OrthogonalMap const& extension =
      message.GetExtension(serialization::OrthogonalMap::extension);
  EXPECT_TRUE(extension.determinant().negative());
  EXPECT_THAT(extension.rotation().quaternion().real_part(),
              AlmostEquals(0.5, 1));
  EXPECT_EQ(0.5,
            extension.rotation().quaternion().imaginary_part().x().double_());
  EXPECT_EQ(0.5,
            extension.rotation().quaternion().imaginary_part().y().double_());
  EXPECT_EQ(0.5,
            extension.rotation().quaternion().imaginary_part().z().double_());
  Orth const o = Orth::ReadFromMessage(message);
  EXPECT_EQ(orthogonal_a_(vector_), o(vector_));
}

}  // namespace internal_orthogonal_map
}  // namespace geometry
}  // namespace principia
