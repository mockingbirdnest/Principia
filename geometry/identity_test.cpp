
#include "geometry/identity.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using quantities::Length;
using quantities::si::Metre;
using ::testing::Eq;

namespace geometry {

class IdentityTest : public testing::Test {
 protected:
  using World1 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST1, Inertial>;
  using World2 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST2, Inertial>;
  using Orth = OrthogonalMap<World1, World2>;
  using Id = Identity<World1, World2>;
  using R3 = R3Element<Length>;

  IdentityTest()
      : vector_(Vector<Length, World1>(
            R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        bivector_(Bivector<Length, World1>(
            R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        trivector_(Trivector<Length, World1>(4.0 * Metre)),
        form_(SymmetricBilinearForm<Length, World1>(
            R3x3Matrix<Length>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre},
                               {2.0 * Metre, 5.0 * Metre, 6.0 * Metre},
                               {3.0 * Metre, 6.0 * Metre, 4.0 * Metre}))) {}

  Vector<Length, World1> const vector_;
  Bivector<Length, World1> const bivector_;
  Trivector<Length, World1> const trivector_;
  SymmetricBilinearForm<Length, World1> const form_;
};

using IdentityDeathTest = IdentityTest;

TEST_F(IdentityTest, Determinant) {
  EXPECT_TRUE(Id().Determinant().is_positive());
}

TEST_F(IdentityTest, AppliedToVector) {
  EXPECT_THAT(Id()(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(IdentityTest, AppliedToBivector) {
  EXPECT_THAT(Id()(bivector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(IdentityTest, AppliedToTrivector) {
  EXPECT_THAT(Id()(trivector_).coordinates(),
              Eq(4.0 * Metre));
}

TEST_F(IdentityTest, AppliedToSymmetricBilinearForm) {
  R3x3Matrix<Length> const expected_coordinates(
      {1.0 * Metre, 2.0 * Metre, 3.0 * Metre},
      {2.0 * Metre, 5.0 * Metre, 6.0 * Metre},
      {3.0 * Metre, 6.0 * Metre, 4.0 * Metre});
  EXPECT_THAT(Id()(form_).coordinates(),
              Eq(expected_coordinates));
}

TEST_F(IdentityTest, Inverse) {
  Vector<Length, World1> const vector1 = vector_;
  Vector<Length, World2> const vector2 =
      Vector<Length, World2>(R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
  EXPECT_THAT(Id().Inverse()(vector2).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
  Id id;
  Identity<World1, World1> const identity1 = id.Inverse() * id;
  EXPECT_THAT(identity1(vector1), Eq(vector1));
  Identity<World2, World2> const identity2 = id * id.Inverse();
  EXPECT_THAT(identity2(vector2), Eq(vector2));
}

TEST_F(IdentityTest, Forget) {
  EXPECT_THAT(Id().Forget()(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(IdentityTest, Compose) {
  using World3 = Frame<enum class World3Tag>;
  using Orth12 = OrthogonalMap<World1, World2>;
  using Orth13 = OrthogonalMap<World1, World3>;
  using Orth23 = OrthogonalMap<World2, World3>;
  using Id12 = Identity<World1, World2>;
  using Id13 = Identity<World1, World3>;
  using Id23 = Identity<World2, World3>;
  Id12 id12;
  Orth12 const o12 = id12.Forget();
  Id23 id23;
  Orth23 const o23 = id23.Forget();
  Id13 const id13 = id23 * id12;
  Orth13 const o13 = o23 * o12;
  for (Length l = 1 * Metre; l < 4 * Metre; l += 1 * Metre) {
    Vector<Length, World1> modified_vector(
        {l, vector_.coordinates().y, vector_.coordinates().z});
    EXPECT_THAT(id13(modified_vector), Eq(o13(modified_vector)));
  }
}

TEST_F(IdentityDeathTest, SerializationError) {
  using Id12 = Identity<World1, World2>;
  EXPECT_DEATH({
    serialization::LinearMap message;
    Id12 const id = Id12::ReadFromMessage(message);
  }, "Fingerprint");
}

TEST_F(IdentityTest, SerializationSuccess) {
  serialization::LinearMap message;
  Identity<World1, World2> id12a;
  id12a.WriteToMessage(&message);
  EXPECT_TRUE(message.has_from_frame());
  EXPECT_TRUE(message.has_to_frame());
  EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
            message.to_frame().tag_type_fingerprint());
  EXPECT_NE(message.from_frame().tag(),
            message.to_frame().tag());
  EXPECT_EQ(message.from_frame().is_inertial(),
            message.to_frame().is_inertial());
  Identity<World1, World2> const id12b =
      Identity<World1, World2>::ReadFromMessage(message);
  EXPECT_THAT(id12a(vector_), id12b(vector_));
}

TEST_F(IdentityTest, Output) {
  using Id12 = Identity<World1, World2>;
  Id12 id12;
  std::cout << id12 << "\n";
}

}  // namespace geometry
}  // namespace principia
