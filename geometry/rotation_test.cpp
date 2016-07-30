
#include "geometry/rotation.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using quantities::si::Degree;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;
using ::testing::Eq;

namespace geometry {

class RotationTest : public testing::Test {
 protected:
  using World =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;
  using World1 =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST1, true>;
  using Orth = OrthogonalMap<World, World>;
  using Rot = Rotation<World, World>;

  RotationTest()
      : vector_(Vector<quantities::Length, World>(
            R3Element<quantities::Length>(
                1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        bivector_(Bivector<quantities::Length, World>(
            R3Element<quantities::Length>(
                1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        trivector_(Trivector<quantities::Length, World>(4.0 * Metre)),
        e1_(Vector<double, World>(R3Element<double>({1, 0, 0}))),
        e2_(Vector<double, World>(R3Element<double>({0, 1, 0}))),
        e3_(Vector<double, World>(R3Element<double>({0, 0, 1}))),
        rotation_a_(Rot(120 * Degree, Bivector<double, World>({1, 1, 1}))),
        rotation_b_(Rot(90 * Degree, Bivector<double, World>({1, 0, 0}))),
        rotation_c_(Rot(ToQuaternion(R3x3Matrix({{0.5, 0.5 * sqrt(3), 0},
                                                 {-0.5 * sqrt(3), 0.5, 0},
                                                 {0, 0, 1}})))) {}

  Vector<quantities::Length, World> vector_;
  Bivector<quantities::Length, World> bivector_;
  Trivector<quantities::Length, World> trivector_;
  Vector<double, World> e1_;
  Vector<double, World> e2_;
  Vector<double, World> e3_;
  Rot rotation_a_;
  Rot rotation_b_;
  Rot rotation_c_;
};

using RotationDeathTest = RotationTest;

TEST_F(RotationTest, Identity) {
  EXPECT_THAT(vector_, Eq(Rot::Identity()(vector_)));
  EXPECT_THAT(bivector_, Eq(Rot::Identity()(bivector_)));
  EXPECT_THAT(trivector_, Eq(Rot::Identity()(trivector_)));
}

TEST_F(RotationTest, AppliedToVector) {
  EXPECT_THAT(rotation_a_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
  EXPECT_THAT(rotation_b_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1, 2));
  EXPECT_THAT(rotation_c_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>((0.5 + sqrt(3.0)) * Metre,
                                                (1.0 - 0.5 * sqrt(3.0)) * Metre,
                                                3.0 * Metre)), 0));
}

TEST_F(RotationTest, AppliedToBivector) {
  EXPECT_THAT(rotation_a_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
  EXPECT_THAT(rotation_b_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1, 2));
  EXPECT_THAT(rotation_c_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>((0.5 + sqrt(3.0)) * Metre,
                                                (1.0 - 0.5 * sqrt(3.0)) * Metre,
                                                3.0 * Metre)), 0));
}

TEST_F(RotationTest, AppliedToTrivector) {
  EXPECT_THAT(rotation_a_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
  EXPECT_THAT(rotation_b_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
  EXPECT_THAT(rotation_c_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
}

TEST_F(RotationTest, Determinant) {
  EXPECT_TRUE(rotation_a_.Determinant().Positive());
  EXPECT_TRUE(rotation_b_.Determinant().Positive());
  EXPECT_TRUE(rotation_c_.Determinant().Positive());
}

TEST_F(RotationTest, Inverse) {
  EXPECT_THAT(rotation_a_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                3.0 * Metre,
                                                1.0 * Metre)), 2));
  EXPECT_THAT(rotation_b_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre)), 1, 2));
  EXPECT_THAT(rotation_c_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>((0.5 - sqrt(3.0)) * Metre,
                                                (1.0 + 0.5 * sqrt(3.0)) * Metre,
                                                3.0 * Metre)), 0));
}

TEST_F(RotationTest, Composition) {
  Rot const rotation_ab = rotation_a_ * rotation_b_;
  EXPECT_THAT(rotation_ab(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                1.0 * Metre,
                                                -3.0 * Metre)), 4, 6));
}

TEST_F(RotationTest, Forget) {
  Orth const orthogonal_a = rotation_a_.Forget();
  EXPECT_THAT(orthogonal_a(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
}

// These four tests cover all the branches of ToQuaternion.
TEST_F(RotationTest, ToQuaternion1) {
  R3Element<double> const v1 = {2, 5, 6};
  R3Element<double> v2 = {-3, 4, 1};
  v1.Orthogonalize<double>(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationTest, ToQuaternion2) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> v2 = {-3, 4, 1};
  v1.Orthogonalize<double>(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationTest, ToQuaternion3) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> v2 = {-3, -4, 1};
  v1.Orthogonalize<double>(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 2));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 1));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 12));
}

TEST_F(RotationTest, ToQuaternion4) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> v2 = {-3, -4, -1};
  v1.Orthogonalize<double>(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 1));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 2));
}

TEST_F(RotationDeathTest, SerializationError) {
  Identity<World, World> id;
  EXPECT_DEATH({
    serialization::LinearMap message;
    id.WriteToMessage(&message);
    Rot const r = Rot::ReadFromMessage(message);
  }, "HasExtension.*Rotation");
}

TEST_F(RotationTest, SerializationSuccess) {
  serialization::LinearMap message;
  rotation_a_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_from_frame());
  EXPECT_TRUE(message.has_to_frame());
  EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
            message.to_frame().tag_type_fingerprint());
  EXPECT_EQ(message.from_frame().tag(),
            message.to_frame().tag());
  EXPECT_EQ(message.from_frame().is_inertial(),
            message.to_frame().is_inertial());
  EXPECT_TRUE(message.HasExtension(serialization::Rotation::extension));
  serialization::Rotation const& extension =
      message.GetExtension(serialization::Rotation::extension);
  EXPECT_THAT(extension.quaternion().real_part(), AlmostEquals(0.5, 1));
  EXPECT_EQ(0.5, extension.quaternion().imaginary_part().x().double_());
  EXPECT_EQ(0.5, extension.quaternion().imaginary_part().y().double_());
  EXPECT_EQ(0.5, extension.quaternion().imaginary_part().z().double_());
  Rot const r = Rot::ReadFromMessage(message);
  EXPECT_EQ(rotation_a_(vector_), r(vector_));
}

TEST_F(RotationTest, Basis) {
  Vector<double, World> a = Normalize(Vector<double, World>({1, 1, -1}));
  Vector<double, World> b = Normalize(Vector<double, World>({1, 0, 1}));
  Bivector<double, World> c = Wedge(a, b);

  Rotation<World, World1> const to_world1(a, b, c);
  EXPECT_THAT(to_world1(a),
              Componentwise(
                  AlmostEquals(1, 1),
                  AlmostEquals(0, 0),
                  VanishesBefore(1, 0)));
  EXPECT_THAT(to_world1(b),
              Componentwise(
                  AlmostEquals(0, 0),
                  AlmostEquals(1, 0),
                  VanishesBefore(1, 0)));
  EXPECT_THAT(to_world1(c),
              Componentwise(
                  VanishesBefore(1, 2),
                  AlmostEquals(0, 0),
                  AlmostEquals(1, 0)));

  Rotation<World1, World> const to_world(a, b, c);
  EXPECT_THAT(to_world(Vector<double, World1>({1, 0, 0})), AlmostEquals(a, 1));
  EXPECT_THAT(to_world(Vector<double, World1>({0, 1, 0})), AlmostEquals(b, 2));
  EXPECT_THAT(to_world(Bivector<double, World1>({0, 0, 1})),
              AlmostEquals(c, 4));
}

}  // namespace geometry
}  // namespace principia
