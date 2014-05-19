#include "geometry/permutation.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

using principia::si::Metre;
using testing::Eq;

namespace principia {
namespace geometry {

class PermutationTest : public testing::Test {
 protected:
  struct World;
  typedef Permutation<World, World> Perm;
  typedef R3Element<quantities::Length> R3;

  void SetUp() override {
    vector_ = Vector<quantities::Length, World>(
        R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    bivector_ = Bivector<quantities::Length, World>(
        R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    trivector_ = Trivector<quantities::Length, World>(4.0 * Metre);
  }

  Vector<quantities::Length, World> vector_;
  Bivector<quantities::Length, World> bivector_;
  Trivector<quantities::Length, World> trivector_;
};

TEST_F(PermutationTest, Identity) {
  EXPECT_THAT(Perm::Identity()(vector_.coordinates()),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, XYZ) {
  EXPECT_THAT(Perm(Perm::XYZ)(vector_.coordinates()),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, YZX) {
  EXPECT_THAT(Perm(Perm::YZX)(vector_.coordinates()),
              Eq<R3>({2.0 * Metre, 3.0 * Metre, 1.0 * Metre}));
}

TEST_F(PermutationTest, ZXY) {
  EXPECT_THAT(Perm(Perm::ZXY)(vector_.coordinates()),
              Eq<R3>({3.0 * Metre, 1.0 * Metre, 2.0 * Metre}));
}

TEST_F(PermutationTest, XZY) {
  EXPECT_THAT(Perm(Perm::XZY)(vector_.coordinates()),
              Eq<R3>({1.0 * Metre, 3.0 * Metre, 2.0 * Metre}));
}

TEST_F(PermutationTest, ZYX) {
  EXPECT_THAT(Perm(Perm::ZYX)(vector_.coordinates()),
              Eq<R3>({3.0 * Metre, 2.0 * Metre, 1.0 * Metre}));
}

TEST_F(PermutationTest, YXZ) {
  EXPECT_THAT(Perm(Perm::YXZ)(vector_.coordinates()),
              Eq<R3>({2.0 * Metre, 1.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, Determinant) {
  Perm xyz(Perm::XYZ);
  Perm yzx(Perm::YZX);
  Perm zxy(Perm::ZXY);
  Perm xzy(Perm::XZY);
  Perm zyx(Perm::ZYX);
  Perm yxz(Perm::YXZ);
  EXPECT_TRUE(xyz.Determinant().Positive());
  EXPECT_TRUE(yzx.Determinant().Positive());
  EXPECT_TRUE(zxy.Determinant().Positive());
  EXPECT_TRUE(xzy.Determinant().Negative());
  EXPECT_TRUE(zyx.Determinant().Negative());
  EXPECT_TRUE(yxz.Determinant().Negative());
}

TEST_F(PermutationTest, AppliedToVector) {
  EXPECT_THAT(Perm(Perm::YZX)(vector_).coordinates(),
              Eq<R3>({2.0 * Metre, 3.0 * Metre, 1.0 * Metre}));
  EXPECT_THAT(Perm(Perm::YXZ)(vector_).coordinates(),
              Eq<R3>({2.0 * Metre, 1.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, AppliedToBivector) {
  EXPECT_THAT(Perm(Perm::ZXY)(bivector_).coordinates(),
              Eq<R3>({3.0 * Metre, 1.0 * Metre, 2.0 * Metre}));
  EXPECT_THAT(Perm(Perm::ZYX)(bivector_).coordinates(),
              Eq<R3>({-3.0 * Metre, -2.0 * Metre, -1.0 * Metre}));
}

TEST_F(PermutationTest, AppliedToTrivector) {
  EXPECT_THAT(Perm(Perm::XYZ)(trivector_).coordinates(),
              Eq(4.0 * Metre));
  EXPECT_THAT(Perm(Perm::XZY)(trivector_).coordinates(),
              Eq(-4.0 * Metre));
}

}  // namespace geometry
}  // namespace principia
