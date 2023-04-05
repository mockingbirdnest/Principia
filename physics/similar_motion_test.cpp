#include "physics/similar_motion.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/homothecy.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/space_transformations.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "geometry/signature.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace physics {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_homothecy;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_space_transformations;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_signature;
using namespace principia::geometry::_space;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_componentwise;

class SimilarMotionTest : public ::testing::Test {
 protected:
  using World1 = Frame<struct World1Tag>;
  using World2 = Frame<struct World2Tag>;
  using World3 = Frame<struct World3Tag>;
  using World4 = Frame<struct World4Tag>;
  using World5 = Frame<struct World5Tag, Arbitrary, Handedness::Left>;
  using World6 = Frame<struct World6Tag, Arbitrary, Handedness::Left>;
  using World7 = Frame<struct World7Tag, Arbitrary, Handedness::Left>;
  using World8 = Frame<struct World8Tag>;
  using World9 = Frame<struct World9Tag>;

  SimilarMotionTest()
      : rotation_a_(
            Vector<double, World1>({Cos(0.2 * Radian), Sin(0.2 * Radian), 0}),
            Vector<double, World1>({-Sin(0.2 * Radian), Cos(0.2 * Radian), 0}),
            Bivector<double, World1>({0, 0, 1})),
        rotation_b_(
            Vector<double, World4>({0, -1, 0}),
            Vector<double, World4>({-Sin(0.8 * Radian), 0, Cos(0.8 * Radian)}),
            Bivector<double, World4>(
                {-Cos(0.8 * Radian), 0, -Sin(0.8 * Radian)})),
        rotation_c_(
            Vector<double, World6>({-1, 0, 0}),
            Vector<double, World6>({0, -Sin(0.7 * Radian), Cos(0.7 * Radian)}),
            Bivector<double, World6>(
                {0, Cos(0.7 * Radian), Sin(0.7 * Radian)})),
        signature_b_(Sign::Negative(),
                     Signature<World4, World5>::DeduceSign(),
                     Sign::Negative()),
        signature_c_(Sign::Positive(),
                     Signature<World7, World8>::DeduceSign(),
                     Sign::Negative()),
        rigid_transformation_a_(
            World1::origin +
                Displacement<World1>({1 * Metre, 0 * Metre, 0 * Metre}),
            World2::origin +
                Displacement<World2>({0 * Metre, 2 * Metre, 0 * Metre}),
            rotation_a_.Forget<OrthogonalMap>()),
        rigid_transformation_b_(
            World3::origin +
                Displacement<World3>({7 * Metre, -8 * Metre, 0 * Metre}),
            World5::origin +
                Displacement<World5>({-1 * Metre, 2 * Metre, 11 * Metre}),
            signature_b_.Forget<OrthogonalMap>() *
                rotation_b_.Forget<OrthogonalMap>()),
        rigid_transformation_c_(
            World6::origin +
                Displacement<World6>({7 * Metre, -3 * Metre, 2 * Metre}),
            World8::origin +
                Displacement<World8>({-13 * Metre, 2 * Metre, 1 * Metre}),
            signature_c_.Forget<OrthogonalMap>() *
                rotation_c_.Forget<OrthogonalMap>()),
        rigid_motion_a_(
            rigid_transformation_a_,
            AngularVelocity<World1>({0 * Radian / Second,
                                     0 * Radian / Second,
                                     3 * Radian / Second}),
            Velocity<World1>(
                {-1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second})),
        rigid_motion_b_(
            rigid_transformation_b_,
            AngularVelocity<World5>({10 * Radian / Second,
                                     3 * Radian / Second,
                                     -3 * Radian / Second}),
            Velocity<World5>(
                {-1 * Metre / Second, 4 * Metre / Second, 0 * Metre / Second})),
        rigid_motion_c_(
            rigid_transformation_c_,
            AngularVelocity<World6>({7 * Radian / Second,
                                     1 * Radian / Second,
                                     -3 * Radian / Second}),
            Velocity<World6>(
                {-5 * Metre / Second, 4 * Metre / Second, 5 * Metre / Second})),
        homothecy_a_(5),
        homothecy_b_(0.25),
        homothecy_c_(0.4),
        similar_motion_a_(SimilarMotion<World2, World3>::DilatationAboutOrigin(
                              homothecy_a_,
                              /*dilatation_rate=*/5 / Second) *
                          rigid_motion_a_.Forget<SimilarMotion>()),
        similar_motion_b_(SimilarMotion<World5, World6>::DilatationAboutOrigin(
                              homothecy_b_,
                              /*dilatation_rate=*/0.3 / Second) *
                          rigid_motion_b_.Forget<SimilarMotion>()),
        similar_motion_c_(SimilarMotion<World8, World9>::DilatationAboutOrigin(
                              homothecy_c_,
                              /*dilatation_rate=*/0.3 / Second) *
                          rigid_motion_c_.Forget<SimilarMotion>()),
        q1_(World1::origin +
            Displacement<World1>({-5 * Metre, 7 * Metre, 11 * Metre})),
        v1_({-2 * Metre / Second, 4 * Metre / Second, -8 * Metre / Second}),
        q6_(World6::origin +
            Displacement<World6>({15 * Metre, 37 * Metre, 101 * Metre})),
        v6_({-2 * Metre / Second, -4 * Metre / Second, -18 * Metre / Second}) {}

  Rotation<World1, World2> const rotation_a_;
  Rotation<World3, World4> const rotation_b_;
  Rotation<World6, World7> const rotation_c_;
  Signature<World4, World5> const signature_b_;
  Signature<World7, World8> const signature_c_;
  RigidTransformation<World1, World2> const rigid_transformation_a_;
  RigidTransformation<World3, World5> const rigid_transformation_b_;
  RigidTransformation<World6, World8> const rigid_transformation_c_;
  RigidMotion<World1, World2> const rigid_motion_a_;
  RigidMotion<World3, World5> const rigid_motion_b_;
  RigidMotion<World6, World8> const rigid_motion_c_;
  Homothecy<double, World2, World3> homothecy_a_;
  Homothecy<double, World5, World6> homothecy_b_;
  Homothecy<double, World8, World9> homothecy_c_;
  SimilarMotion<World1, World3> const similar_motion_a_;
  SimilarMotion<World3, World6> const similar_motion_b_;
  SimilarMotion<World6, World9> const similar_motion_c_;

  Position<World1> q1_;
  Velocity<World1> v1_;
  Position<World6> q6_;
  Velocity<World6> v6_;
};

TEST_F(SimilarMotionTest, Inverse) {
  EXPECT_THAT((similar_motion_a_.Inverse() * similar_motion_a_)({q1_, v1_}),
              Componentwise(AlmostEquals(q1_, 0), AlmostEquals(v1_, 5)));
  EXPECT_THAT(similar_motion_a_.Inverse()(similar_motion_a_({q1_, v1_})),
              Componentwise(AlmostEquals(q1_, 0), AlmostEquals(v1_, 24)));
  EXPECT_THAT((similar_motion_b_ * similar_motion_b_.Inverse())({q6_, v6_}),
              Componentwise(AlmostEquals(q6_, 0), AlmostEquals(v6_, 3)));
  EXPECT_THAT(similar_motion_b_(similar_motion_b_.Inverse()({q6_, v6_})),
              Componentwise(AlmostEquals(q6_, 8), AlmostEquals(v6_, 512)));
}

TEST_F(SimilarMotionTest, Composition) {
  auto const qv6 = similar_motion_b_(similar_motion_a_({q1_, v1_}));
  EXPECT_THAT((similar_motion_b_ * similar_motion_a_)({q1_, v1_}),
              Componentwise(AlmostEquals(qv6.position(), 7),
                            AlmostEquals(qv6.velocity(), 3)));
}

TEST_F(SimilarMotionTest, Associativity) {
  auto const qv9 =
      similar_motion_c_(similar_motion_b_(similar_motion_a_({q1_, v1_})));
  EXPECT_THAT(
      ((similar_motion_c_ * similar_motion_b_) * similar_motion_a_)({q1_, v1_}),
      Componentwise(AlmostEquals(qv9.position(), 8),
                    AlmostEquals(qv9.velocity(), 12)));
  EXPECT_THAT(
      (similar_motion_c_ * (similar_motion_b_ * similar_motion_a_))({q1_, v1_}),
      Componentwise(AlmostEquals(qv9.position(), 4),
                    AlmostEquals(qv9.velocity(), 6)));
  EXPECT_THAT(
      (similar_motion_c_ * similar_motion_b_)(similar_motion_a_({ q1_, v1_ })),
      Componentwise(AlmostEquals(qv9.position(), 4),
                    AlmostEquals(qv9.velocity(), 7)));
  EXPECT_THAT(
      similar_motion_c_((similar_motion_b_ * similar_motion_a_)({q1_, v1_})),
      Componentwise(AlmostEquals(qv9.position(), 5),
                    AlmostEquals(qv9.velocity(), 5)));
}

TEST_F(SimilarMotionTest, Forget) {
  auto const qv2 = rigid_motion_a_({q1_, v1_});
  EXPECT_THAT(rigid_motion_a_.Forget<SimilarMotion>()({q1_, v1_}),
              Componentwise(AlmostEquals(qv2.position(), 0),
                            AlmostEquals(qv2.velocity(), 1)));
}

}  // namespace physics
}  // namespace principia
