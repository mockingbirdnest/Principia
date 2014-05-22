#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using si::Metre;
using testing::Eq;
using testing_utilities::AlmostEquals;

class RotationTest : public testing::Test {
 protected:
  struct World;
  typedef OrthogonalMap<World, World> Orth;
  typedef Rotation<World, World> Rot;

  void SetUp() override {
    vector_ = Vector<quantities::Length, World>(
        R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    bivector_ = Bivector<quantities::Length, World>(
        R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    trivector_ = Trivector<quantities::Length, World>(4.0 * Metre);
    rotation_a_ = Rot(120 * si::Degree,
                      Vector<quantities::Dimensionless, World>({1, 1, 1}));
    rotation_b_ = Rot(90 * si::Degree,
                      Vector<quantities::Dimensionless, World>({1, 0, 0}));
  }

  Vector<quantities::Length, World> vector_;
  Bivector<quantities::Length, World> bivector_;
  Trivector<quantities::Length, World> trivector_;
  Rot rotation_a_;
  Rot rotation_b_;
};

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
                                                2.0 * Metre))));
  EXPECT_THAT(rotation_b_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre))));
}

TEST_F(RotationTest, AppliedToBivector) {
  EXPECT_THAT(rotation_a_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre))));
  EXPECT_THAT(rotation_b_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre))));
}

TEST_F(RotationTest, AppliedToTrivector) {
  EXPECT_THAT(rotation_a_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre)));
  EXPECT_THAT(rotation_b_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre)));
}

TEST_F(RotationTest, Determinant) {
  EXPECT_TRUE(rotation_a_.Determinant().Positive());
  EXPECT_TRUE(rotation_b_.Determinant().Positive());
}

TEST_F(RotationTest, Inverse) {
  EXPECT_THAT(rotation_a_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                3.0 * Metre,
                                                1.0 * Metre))));
  EXPECT_THAT(rotation_b_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre))));
}

TEST_F(RotationTest, Composition) {
  Rot const rotation_ab = rotation_a_ * rotation_b_;
  EXPECT_THAT(rotation_ab(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                1.0 * Metre,
                                                -3.0 * Metre))));
}

TEST_F(RotationTest, Forget) {
  Orth const orthogonal_a = rotation_a_.Forget();
  EXPECT_THAT(rotation_a_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre))));
}

}  // namespace geometry
}  // namespace principia
