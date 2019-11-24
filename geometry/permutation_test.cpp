
#include "geometry/permutation.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/identity.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using quantities::Length;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;
using ::testing::Eq;

namespace geometry {

class PermutationTest : public testing::Test {
 protected:
  using World1 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST1, true>;
  using World2 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST2, true>;
  using Orth = OrthogonalMap<World1, World2>;
  using Perm = Permutation<World1, World2>;
  using R3 = R3Element<quantities::Length>;

  PermutationTest() {
    vector_ = Vector<quantities::Length, World1>(
        R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    bivector_ = Bivector<quantities::Length, World1>(
        R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    trivector_ = Trivector<quantities::Length, World1>(4.0 * Metre);
  }

  Vector<quantities::Length, World1> vector_;
  Bivector<quantities::Length, World1> bivector_;
  Trivector<quantities::Length, World1> trivector_;
};

using PermutationDeathTest = PermutationTest;

TEST_F(PermutationTest, Identity) {
  EXPECT_THAT(Perm::Identity()(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, XYZ) {
  EXPECT_THAT(Perm(Perm::XYZ)(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, YZX) {
  EXPECT_THAT(Perm(Perm::YZX)(vector_).coordinates(),
              Eq<R3>({2.0 * Metre, 3.0 * Metre, 1.0 * Metre}));
}

TEST_F(PermutationTest, ZXY) {
  EXPECT_THAT(Perm(Perm::ZXY)(vector_).coordinates(),
              Eq<R3>({3.0 * Metre, 1.0 * Metre, 2.0 * Metre}));
}

TEST_F(PermutationTest, XZY) {
  EXPECT_THAT(Perm(Perm::XZY)(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 3.0 * Metre, 2.0 * Metre}));
}

TEST_F(PermutationTest, ZYX) {
  EXPECT_THAT(Perm(Perm::ZYX)(vector_).coordinates(),
              Eq<R3>({3.0 * Metre, 2.0 * Metre, 1.0 * Metre}));
}

TEST_F(PermutationTest, YXZ) {
  EXPECT_THAT(Perm(Perm::YXZ)(vector_).coordinates(),
              Eq<R3>({2.0 * Metre, 1.0 * Metre, 3.0 * Metre}));
}

TEST_F(PermutationTest, Determinant) {
  Perm xyz(Perm::XYZ);
  Perm yzx(Perm::YZX);
  Perm zxy(Perm::ZXY);
  Perm xzy(Perm::XZY);
  Perm zyx(Perm::ZYX);
  Perm yxz(Perm::YXZ);
  EXPECT_TRUE(xyz.Determinant().is_positive());
  EXPECT_TRUE(yzx.Determinant().is_positive());
  EXPECT_TRUE(zxy.Determinant().is_positive());
  EXPECT_TRUE(xzy.Determinant().is_negative());
  EXPECT_TRUE(zyx.Determinant().is_negative());
  EXPECT_TRUE(yxz.Determinant().is_negative());
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
  Vector<quantities::Length, World1> const vector1 = vector_;
  Vector<quantities::Length, World2> const vector2 =
      Vector<quantities::Length, World2>(
          R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
  EXPECT_THAT(Perm(Perm::YZX).Inverse()(vector2).coordinates(),
              Eq<R3>({3.0 * Metre, 1.0 * Metre, 2.0 * Metre}));
  EXPECT_THAT(Perm(Perm::YXZ).Inverse()(vector2).coordinates(),
              Eq<R3>({2.0 * Metre, 1.0 * Metre, 3.0 * Metre}));

  std::vector<Perm> const all =
      {Perm(Perm::XYZ), Perm(Perm::YZX), Perm(Perm::ZXY),
       Perm(Perm::XZY), Perm(Perm::ZYX), Perm(Perm::YXZ)};
  for (Perm const& p : all) {
    Permutation<World1, World1> const identity1 = p.Inverse() * p;
    EXPECT_THAT(identity1(vector1), Eq(vector1));
    Permutation<World2, World2> const identity2 = p * p.Inverse();
    EXPECT_THAT(identity2(vector2), Eq(vector2));
  }
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
  struct World3;
  using Orth12 = OrthogonalMap<World1, World2>;
  using Orth13 = OrthogonalMap<World1, World3>;
  using Orth23 = OrthogonalMap<World2, World3>;
  using Perm12 = Permutation<World1, World2>;
  using Perm13 = Permutation<World1, World3>;
  using Perm23 = Permutation<World2, World3>;
  std::vector<Perm12> const all12 =
      {Perm12(Perm12::XYZ), Perm12(Perm12::YZX), Perm12(Perm12::ZXY),
       Perm12(Perm12::XZY), Perm12(Perm12::ZYX), Perm12(Perm12::YXZ)};
  std::vector<Perm23> const all23 =
      {Perm23(Perm23::XYZ), Perm23(Perm23::YZX), Perm23(Perm23::ZXY),
       Perm23(Perm23::XZY), Perm23(Perm23::ZYX), Perm23(Perm23::YXZ)};
  for (Perm12 const& p12 : all12) {
    Orth12 const o12 = p12.Forget();
    for (Perm23 const& p23 : all23) {
      Orth23 const o23 = p23.Forget();
      Perm13 const p13 = p23 * p12;
      Orth13 const o13 = o23 * o12;
      for (Length l = 1 * Metre; l < 4 * Metre; l += 1 * Metre) {
        Vector<quantities::Length, World1> modified_vector(
            {l, vector_.coordinates().y, vector_.coordinates().z});
        EXPECT_THAT(p13(modified_vector),
                    AlmostEquals(o13(modified_vector), 0, 12));
      }
    }
  }
}

TEST_F(PermutationDeathTest, SerializationError) {
  Identity<World1, World2> id;
  EXPECT_DEATH({
    serialization::LinearMap message;
    id.WriteToMessage(&message);
    Perm const p = Perm::ReadFromMessage(message);
  }, "HasExtension.*Permutation");
}

TEST_F(PermutationTest, SerializationSuccess) {
  std::vector<Perm::CoordinatePermutation> const all12 =
      {Perm::XYZ, Perm::YZX, Perm::ZXY,
       Perm::XZY, Perm::ZYX, Perm::YXZ};
  serialization::LinearMap message;
  for (Perm::CoordinatePermutation const cp : all12) {
    Perm const perm_a(cp);
    perm_a.WriteToMessage(&message);
    EXPECT_TRUE(message.has_from_frame());
    EXPECT_TRUE(message.has_to_frame());
    EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
              message.to_frame().tag_type_fingerprint());
    EXPECT_NE(message.from_frame().tag(),
              message.to_frame().tag());
    EXPECT_EQ(message.from_frame().is_inertial(),
              message.to_frame().is_inertial());
    EXPECT_TRUE(message.HasExtension(serialization::Permutation::extension));
    serialization::Permutation const& extension =
        message.GetExtension(serialization::Permutation::extension);
    EXPECT_EQ(cp, extension.coordinate_permutation());
    Perm const perm_b = Perm::ReadFromMessage(message);
    EXPECT_EQ(perm_a(vector_), perm_b(vector_));
  }
}

}  // namespace geometry
}  // namespace principia
