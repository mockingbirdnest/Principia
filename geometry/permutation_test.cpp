
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
#include "testing_utilities/componentwise.hpp"

namespace principia {

using quantities::Length;
using quantities::si::Metre;
using ::testing::Eq;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;

namespace geometry {

class PermutationTest : public testing::Test {
 protected:
  using R1 = Frame<serialization::Frame::TestTag,
                   Inertial,
                   Handedness::Right,
                   serialization::Frame::TEST1>;
  using R2 = Frame<serialization::Frame::TestTag,
                   Inertial,
                   Handedness::Right,
                   serialization::Frame::TEST2>;
  using L = Frame<serialization::Frame::TestTag,
                  Inertial,
                  Handedness::Left,
                  serialization::Frame::TEST3>;
  using PermutationR1R2 = Permutation<R1, R2>;
  using PermutationR1L = Permutation<R1, L>;

  PermutationTest()
      : vector_({1 * Metre, 2 * Metre, 3 * Metre}),
        bivector_({1 * Metre, 2 * Metre, 3 * Metre}),
        trivector_(4 * Metre) {}

  Vector<Length, R1> const vector_;
  Bivector<Length, R1> const bivector_;
  Trivector<Length, R1> const trivector_;
};

using PermutationDeathTest = PermutationTest;

TEST_F(PermutationTest, Identity) {
  EXPECT_THAT(PermutationR1R2::Identity()(vector_),
              Componentwise(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
}

TEST_F(PermutationTest, XYZ) {
  EXPECT_THAT(PermutationR1R2(EvenPermutation::XYZ)(vector_),
              Componentwise(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
}

TEST_F(PermutationTest, YZX) {
  EXPECT_THAT(PermutationR1R2(EvenPermutation::YZX)(vector_),
              Componentwise(2.0 * Metre, 3.0 * Metre, 1.0 * Metre));
}

TEST_F(PermutationTest, ZXY) {
  EXPECT_THAT(PermutationR1R2(EvenPermutation::ZXY)(vector_),
              Componentwise(3.0 * Metre, 1.0 * Metre, 2.0 * Metre));
}

TEST_F(PermutationTest, XZY) {
  EXPECT_THAT(PermutationR1L(OddPermutation::XZY)(vector_),
              Componentwise(1.0 * Metre, 3.0 * Metre, 2.0 * Metre));
}

TEST_F(PermutationTest, ZYX) {
  EXPECT_THAT(PermutationR1L(OddPermutation::ZYX)(vector_),
              Componentwise(3.0 * Metre, 2.0 * Metre, 1.0 * Metre));
}

TEST_F(PermutationTest, YXZ) {
  EXPECT_THAT(PermutationR1L(OddPermutation::YXZ)(vector_),
              Componentwise(2.0 * Metre, 1.0 * Metre, 3.0 * Metre));
}

TEST_F(PermutationTest, Determinant) {
  PermutationR1R2 xyz(EvenPermutation::XYZ);
  PermutationR1R2 yzx(EvenPermutation::YZX);
  PermutationR1R2 zxy(EvenPermutation::ZXY);
  PermutationR1L xzy(OddPermutation::XZY);
  PermutationR1L zyx(OddPermutation::ZYX);
  PermutationR1L yxz(OddPermutation::YXZ);
  EXPECT_TRUE(xyz.Determinant().is_positive());
  EXPECT_TRUE(yzx.Determinant().is_positive());
  EXPECT_TRUE(zxy.Determinant().is_positive());
  EXPECT_TRUE(xzy.Determinant().is_negative());
  EXPECT_TRUE(zyx.Determinant().is_negative());
  EXPECT_TRUE(yxz.Determinant().is_negative());
}

TEST_F(PermutationTest, AppliedToVector) {
  EXPECT_THAT(PermutationR1R2(EvenPermutation::YZX)(vector_),
              Componentwise(2.0 * Metre, 3.0 * Metre, 1.0 * Metre));
  EXPECT_THAT(PermutationR1L(OddPermutation::YXZ)(vector_),
              Componentwise(2.0 * Metre, 1.0 * Metre, 3.0 * Metre));
}

TEST_F(PermutationTest, AppliedToBivector) {
  EXPECT_THAT(PermutationR1R2(EvenPermutation::ZXY)(bivector_),
              Componentwise(3.0 * Metre, 1.0 * Metre, 2.0 * Metre));
  EXPECT_THAT(PermutationR1L(OddPermutation::ZYX)(bivector_),
              Componentwise(-3.0 * Metre, -2.0 * Metre, -1.0 * Metre));
}

TEST_F(PermutationTest, AppliedToTrivector) {
  EXPECT_THAT(PermutationR1R2(EvenPermutation::XYZ)(trivector_).coordinates(),
              Eq(4.0 * Metre));
  EXPECT_THAT(PermutationR1L(OddPermutation::XZY)(trivector_).coordinates(),
              Eq(-4.0 * Metre));
}

TEST_F(PermutationTest, Inverse) {
  EXPECT_THAT(PermutationR1R2(EvenPermutation::YZX)
                  .Inverse()(Vector<Length, R2>(vector_.coordinates())),
              Componentwise(3.0 * Metre, 1.0 * Metre, 2.0 * Metre));
  EXPECT_THAT(PermutationR1L(OddPermutation::YXZ)
                  .Inverse()(Vector<Length, L>(vector_.coordinates())),
              Componentwise(2.0 * Metre, 1.0 * Metre, 3.0 * Metre));

  std::array<PermutationR1R2, 3> const all_even_permutations = {
      PermutationR1R2(EvenPermutation::XYZ),
      PermutationR1R2(EvenPermutation::YZX),
      PermutationR1R2(EvenPermutation::ZXY)};
  std::array<PermutationR1L, 3> const all_odd_permutations = {
      PermutationR1L(OddPermutation::XZY),
      PermutationR1L(OddPermutation::ZYX),
      PermutationR1L(OddPermutation::YXZ)};
  auto const test_inverse_for = [&](auto const& permutations) {
    for (auto const& p : permutations) {
      auto const identity1 = p.Inverse() * p;
      EXPECT_THAT(identity1(vector_), Eq(vector_));
      auto const identity2 = p * p.Inverse();
      auto const vector_to = decltype(p(vector_))(vector_.coordinates());
      EXPECT_THAT(identity2(vector_to), Eq(vector_to));
    }
  };
  test_inverse_for(all_even_permutations);
  test_inverse_for(all_odd_permutations);
}

TEST_F(PermutationTest, Forget) {
  EXPECT_THAT(
      PermutationR1R2(EvenPermutation::XYZ).Forget<OrthogonalMap>()(vector_),
      Componentwise(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
  EXPECT_THAT(
      PermutationR1R2(EvenPermutation::YZX).Forget<OrthogonalMap>()(vector_),
      Componentwise(2.0 * Metre, 3.0 * Metre, 1.0 * Metre));
  EXPECT_THAT(
      PermutationR1R2(EvenPermutation::ZXY).Forget<OrthogonalMap>()(vector_),
      Componentwise(3.0 * Metre, 1.0 * Metre, 2.0 * Metre));
  EXPECT_THAT(
      PermutationR1L(OddPermutation::XZY).Forget<OrthogonalMap>()(vector_),
      Componentwise(AlmostEquals(1.0 * Metre, 2),
                    AlmostEquals(3.0 * Metre, 2),
                    AlmostEquals(2.0 * Metre, 2)));
  EXPECT_THAT(
      PermutationR1L(OddPermutation::ZYX).Forget<OrthogonalMap>()(vector_),
      Componentwise(AlmostEquals(3.0 * Metre, 2),
                    AlmostEquals(2.0 * Metre, 2),
                    AlmostEquals(1.0 * Metre, 4)));
  EXPECT_THAT(
      PermutationR1L(OddPermutation::YXZ).Forget<OrthogonalMap>()(vector_),
      Componentwise(AlmostEquals(2.0 * Metre, 1),
                    AlmostEquals(1.0 * Metre, 2),
                    AlmostEquals(3.0 * Metre, 2)));
}

TEST_F(PermutationTest, Compose) {
  using Orth12 = OrthogonalMap<R1, R2>;
  using Orth1L = OrthogonalMap<R1, L>;
  using Orth2L = OrthogonalMap<R2, L>;
  using Perm12 = Permutation<R1, R2>;
  using PermL2 = Permutation<L, R2>;
  using Perm1L = Permutation<R1, L>;
  using Perm2L = Permutation<R2, L>;
  using Perm21 = Permutation<R2, R1>;
  std::array<Perm12, 3> const all_12 = {Perm12(EvenPermutation::XYZ),
                                        Perm12(EvenPermutation::YZX),
                                        Perm12(EvenPermutation::ZXY)};
  std::array<PermL2, 3> const all_l2 = {PermL2(OddPermutation::XZY),
                                        PermL2(OddPermutation::ZYX),
                                        PermL2(OddPermutation::YXZ)};
  std::array<Perm21, 3> const all_21 = {Perm21(EvenPermutation::XYZ),
                                        Perm21(EvenPermutation::YZX),
                                        Perm21(EvenPermutation::ZXY)};
  std::array<Perm2L, 3> const all_2l = {Perm2L(OddPermutation::XZY),
                                        Perm2L(OddPermutation::ZYX),
                                        Perm2L(OddPermutation::YXZ)};
  auto const test_composition_for = [&](auto const& lhs,
                                        auto const& rhs,
                                        auto const from_vector) {
    for (auto const& left : lhs) {
      for (auto const& right : rhs) {
        auto const composition = left * right;
        auto const composition_as_orthogonal_maps =
            left.Forget<OrthogonalMap>() * right.Forget<OrthogonalMap>();
        for (Length l = 1 * Metre; l < 4 * Metre; l += 1 * Metre) {
          // TODO(egg): In C++20 we could have template parameters on the lambda
          // which would allow us to deduce this type, instead of having to pass
          // an otherwise unused value which we feed to the decltype.
          decltype(from_vector) const modified_vector(
              {l, vector_.coordinates().y, vector_.coordinates().z});
          EXPECT_THAT(
              composition(modified_vector),
              AlmostEquals(
                  composition_as_orthogonal_maps(modified_vector), 0, 12));
        }
      }
    }
  };
  test_composition_for(all_21, all_12, Vector<Length, R1>{});
  test_composition_for(all_2l, all_12, Vector<Length, R1>{});
  test_composition_for(all_21, all_l2, Vector<Length, L>{});
  test_composition_for(all_2l, all_l2, Vector<Length, L>{});
}

TEST_F(PermutationDeathTest, SerializationError) {
  Identity<R1, R2> id;
  EXPECT_DEATH({
    serialization::LinearMap message;
    id.WriteToMessage(&message);
    PermutationR1R2 const p = PermutationR1R2::ReadFromMessage(message);
  }, "HasExtension.*Permutation");
}

TEST_F(PermutationTest, SerializationSuccess) {
  std::array<EvenPermutation, 3> const all_even_permutations = {
      EvenPermutation::XYZ, EvenPermutation::YZX, EvenPermutation::ZXY};
  std::array<OddPermutation, 3> const all_odd_permutations = {
      OddPermutation::XZY, OddPermutation::ZYX, OddPermutation::YXZ};
  auto const test_serialization_for = [&](auto const& permutations) {
    serialization::LinearMap message;
    for (auto const cp : permutations) {
      using Perm = std::conditional_t<
          std::is_same_v<decltype(cp), EvenPermutation const>,
                             PermutationR1R2,
                             PermutationR1L>;
      Perm const perm_a(cp);
      perm_a.WriteToMessage(&message);
      EXPECT_TRUE(message.has_from_frame());
      EXPECT_TRUE(message.has_to_frame());
      EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
                message.to_frame().tag_type_fingerprint());
      EXPECT_NE(message.from_frame().tag(), message.to_frame().tag());
      EXPECT_EQ(message.from_frame().is_inertial(),
                message.to_frame().is_inertial());
      EXPECT_TRUE(message.HasExtension(serialization::Permutation::extension));
      serialization::Permutation const& extension =
          message.GetExtension(serialization::Permutation::extension);
      EXPECT_EQ(static_cast<int>(cp), extension.coordinate_permutation());
      Perm const perm_b = Perm::ReadFromMessage(message);
      EXPECT_EQ(perm_a(vector_), perm_b(vector_));
    }
  };
  test_serialization_for(all_even_permutations);
  test_serialization_for(all_odd_permutations);
}

TEST_F(PermutationTest, Output) {
  std::array<PermutationR1R2, 3> const all_even_permutations = {
      PermutationR1R2(EvenPermutation::XYZ),
      PermutationR1R2(EvenPermutation::YZX),
      PermutationR1R2(EvenPermutation::ZXY)};
  std::array<PermutationR1L, 3> const all_odd_permutations = {
      PermutationR1L(OddPermutation::XZY),
      PermutationR1L(OddPermutation::ZYX),
      PermutationR1L(OddPermutation::YXZ)};
  auto const test_output_for = [&](auto const& permutations) {
    for (auto const& p : permutations) {
      std::cout << p;
    }
  };
  test_output_for(all_even_permutations);
  test_output_for(all_odd_permutations);
}

}  // namespace geometry
}  // namespace principia
