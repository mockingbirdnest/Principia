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

using namespace si;
using namespace testing;
using namespace testing_utilities;

class OrthogonalMapTest : public testing::Test {
 protected:
  struct World;
  typedef OrthogonalMap<World, World> O;
  typedef Rotation<World, World> R;

  void SetUp() override {
    vector_ = Vector<quantities::Length, World>(
        R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    bivector_ = Bivector<quantities::Length, World>(
        R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    trivector_ = Trivector<quantities::Length, World>(4.0 * Metre);
    orthogonal_a_ = O(Sign(-1),
                      R(120 * si::Degree, 
                        Vector<quantities::Dimensionless, World>({1, 1, 1})));
    orthogonal_b_ = O(Sign(1),
                      R(90 * si::Degree,
                        Vector<quantities::Dimensionless, World>({1, 0, 0})));
    orthogonal_c_ = O(Sign(-1),
                      R(90 * si::Degree,
                        Vector<quantities::Dimensionless, World>({1, 0, 0})));
  }

  Vector<quantities::Length, World> vector_;
  Bivector<quantities::Length, World> bivector_;
  Trivector<quantities::Length, World> trivector_;
  O orthogonal_a_;
  O orthogonal_b_;
  O orthogonal_c_;
};

TEST_F(OrthogonalMapTest, Identity) {
  EXPECT_THAT(vector_, Eq(O::Identity()(vector_)));
  EXPECT_THAT(bivector_, Eq(O::Identity()(bivector_)));
  EXPECT_THAT(trivector_, Eq(O::Identity()(trivector_)));
}

TEST_F(OrthogonalMapTest, AppliedToVector) {
  EXPECT_THAT(orthogonal_a_(vector_), 
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(-3.0 * Metre,
                                                -1.0 * Metre,
                                                -2.0 * Metre))));
  EXPECT_THAT(orthogonal_b_(vector_), 
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre))));
}

TEST_F(OrthogonalMapTest, AppliedToBivector) {
  EXPECT_THAT(orthogonal_a_(bivector_), 
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre))));
  EXPECT_THAT(orthogonal_b_(vector_), 
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre))));
}

TEST_F(OrthogonalMapTest, AppliedToTrivector) {
  EXPECT_THAT(orthogonal_a_(trivector_), 
              AlmostEquals(Trivector<quantities::Length, World>(
                  -4.0 * Metre)));
  EXPECT_THAT(orthogonal_b_(trivector_), 
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre)));
}

TEST_F(OrthogonalMapTest, Determinant) {
  EXPECT_TRUE(orthogonal_a_.Determinant().Negative());
  EXPECT_TRUE(orthogonal_b_.Determinant().Positive());
  EXPECT_TRUE(orthogonal_c_.Determinant().Negative());
}

TEST_F(OrthogonalMapTest, Inverse) {
  EXPECT_THAT(orthogonal_a_.Inverse()(vector_), 
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(-2.0 * Metre,
                                                -3.0 * Metre,
                                                -1.0 * Metre))));
  EXPECT_THAT(orthogonal_b_.Inverse()(vector_), 
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre))));
}

TEST_F(OrthogonalMapTest, Composition) {
  O const orthogonal_ac = orthogonal_a_ * orthogonal_c_;
  EXPECT_THAT(orthogonal_ac(vector_), 
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                1.0 * Metre,
                                                -3.0 * Metre))));
}

}  // namespace geometry
}  // namespace principia
