#include "geometry/identity.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

using principia::quantities::Length;
using principia::si::Metre;
using testing::Eq;

namespace principia {
namespace geometry {

class IdentityTest : public testing::Test {
 protected:
  using World1 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST1, true>;
  using World2 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST2, true>;
  using Orth = OrthogonalMap<World1, World2>;
  using Id = Identity<World1, World2>;
  using R3 = R3Element<quantities::Length>;

  void SetUp() override {
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

using IdentityDeathTest = IdentityTest;

TEST_F(IdentityTest, Determinant) {
  Id identity;
  EXPECT_TRUE(identity.Determinant().Positive());
}

TEST_F(IdentityTest, AppliedToVector) {
  EXPECT_THAT(Id::Identity()(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(IdentityTest, AppliedToBivector) {
  EXPECT_THAT(Id::Identity()(bivector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(IdentityTest, AppliedToTrivector) {
  EXPECT_THAT(Id::Identity()(trivector_).coordinates(),
              Eq(4.0 * Metre));
}

TEST_F(IdentityTest, Inverse) {
  Vector<quantities::Length, World1> const vector1 = vector_;
  Vector<quantities::Length, World2> const vector2 =
      Vector<quantities::Length, World2>(
          R3(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
  EXPECT_THAT(Id::Identity().Inverse()(vector2).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
  Id id;
  Identity<World1, World1> const identity1 = id.Inverse() * id;
  EXPECT_THAT(identity1(vector1), Eq(vector1));
  Identity<World2, World2> const identity2 = id * id.Inverse();
  EXPECT_THAT(identity2(vector2), Eq(vector2));
}

TEST_F(IdentityTest, Forget) {
  EXPECT_THAT(Id::Identity().Forget()(vector_).coordinates(),
              Eq<R3>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre}));
}

TEST_F(IdentityTest, Compose) {
  struct World3;
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
    Vector<quantities::Length, World1> modified_vector(
        {l, vector_.coordinates().y, vector_.coordinates().z});
    EXPECT_THAT(id13(modified_vector), Eq(o13(modified_vector)));
  }
}

TEST_F(IdentityDeathTest, SerializationError) {
  using Id12 = Identity<World1, World2>;
  EXPECT_DEATH({
    serialization::LinearMap message;
    Id12 const id = Id12::ReadFromMessage(message);
  }, "HasExtension.*Identity");
}

TEST_F(IdentityTest, SerializationSuccess) {
  serialization::LinearMap message;
  Identity<World1, World2> id12a;
  id12a.WriteToMessage(&message);
  Identity<World1, World2> const id12b =
      Identity<World1, World2>::ReadFromMessage(message);
  EXPECT_THAT(id12a(vector_), id12b(vector_));
}

}  // namespace geometry
}  // namespace principia
