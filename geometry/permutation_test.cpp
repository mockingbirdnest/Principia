
#include <vector>

#include "geometry/permutation.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

using principia::si::Metre;
using principia::testing_utilities::AlmostEquals;
using testing::Eq;

namespace principia {
namespace geometry {

class PermutationTest : public testing::Test {
 protected:
  struct World;
  typedef OrthogonalMap<World, World> Orth;
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

TEST_F(PermutationTest, Inverse) {
  EXPECT_THAT(Perm(Perm::YZX).Inverse()(vector_).coordinates(),
              Eq<R3>({3.0 * Metre, 1.0 * Metre, 2.0 * Metre}));
  EXPECT_THAT(Perm(Perm::YXZ).Inverse()(vector_).coordinates(),
              Eq<R3>({2.0 * Metre, 1.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, Forget) {
  EXPECT_THAT(Perm(Perm::XYZ).Forget()(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
  EXPECT_THAT(Perm(Perm::YZX).Forget()(vector_).coordinates(),
              Eq<R3>({2.0 * Metre, 3.0 * Metre, 1.0 * Metre}));
  EXPECT_THAT(Perm(Perm::ZXY).Forget()(vector_).coordinates(),
              Eq<R3>({3.0 * Metre, 1.0 * Metre, 2.0 * Metre}));
  EXPECT_THAT(Perm(Perm::XZY).Forget()(vector_).coordinates(),
              AlmostEquals<R3>({1.0 * Metre, 3.0 * Metre, 2.0 * Metre}, 2));
  EXPECT_THAT(Perm(Perm::ZYX).Forget()(vector_).coordinates(),
              AlmostEquals<R3>({3.0 * Metre, 2.0 * Metre, 1.0 * Metre}, 4));
  EXPECT_THAT(Perm(Perm::YXZ).Forget()(vector_).coordinates(),
              AlmostEquals<R3>({2.0 * Metre, 1.0 * Metre, 3.0 * Metre}, 2));
}

TEST_F(PermutationTest, Compose) {
  std::vector<Perm> all = {Perm(Perm::XYZ),Perm(Perm::YZX),Perm(Perm::ZXY),
                           Perm(Perm::XZY),Perm(Perm::ZYX),Perm(Perm::YXZ)};
  for (const Perm& p1 : all) {
    Orth const o1 = p1.Forget();
    for (const Perm& p2 : all) {
      Orth const o2 = p2.Forget();
      Perm const p12 = p1 * p2;
      Orth const o12 = o1 * o2;
      EXPECT_THAT(p12(vector_), AlmostEquals(o12(vector_)));
    }
  }
}

}  // namespace geometry
}  // namespace principia
