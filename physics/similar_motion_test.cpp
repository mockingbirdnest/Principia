#include "physics/similar_motion.hpp"

#include "geometry/frame.hpp"
#include "geometry/homothecy.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rigid_transformation.hpp"
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
using namespace principia::geometry::_homothecy;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_rigid_transformation;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_signature;
using namespace principia::geometry::_space;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_similar_motion;
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
};

TEST_F(SimilarMotionTest, Smoke) {
  Signature<World1, World2> const signature(
      Sign::Negative(),
      Signature<World1, World1>::DeduceSign(),
      Sign::Negative());
  RigidTransformation<World1, World2> const rigid_transformation(
      World1::origin + Displacement<World1>({1 * Metre, 0 * Metre, 0 * Metre}),
      World2::origin + Displacement<World2>({0 * Metre, 2 * Metre, 0 * Metre}),
      signature.Forget<OrthogonalMap>());
  RigidMotion<World1, World2> const rigid_motion(
      rigid_transformation,
      AngularVelocity<World1>(
          {0 * Radian / Second, 0 * Radian / Second, 3 * Radian / Second}),
      Velocity<World1>(
          {-1 * Metre / Second, 0 * Metre / Second, 0 * Metre / Second}));
  Homothecy<double, World2, World3> homothecy =
      Homothecy<Length, World4, World3>(4 * Metre) *
      Homothecy<Inverse<Length>, World2, World4>(1 / Metre);
  SimilarMotion<World1, World3> const similar_motion =
      SimilarMotion<World2, World3>::DilatationAboutOrigin(
          homothecy, /*dilatation_rate=*/5 / Second) *
      rigid_motion.Forget<SimilarMotion>();

  auto const transformed_unmoving_origin =
      similar_motion({World1::origin, World1::unmoving});
  EXPECT_THAT(
      transformed_unmoving_origin,
      Componentwise(
          AlmostEquals(
              World3::origin + Displacement<World3>({4 * Metre,
                                                     8 * Metre,
                                                     0 * Metre}), 0),
          AlmostEquals(
              World3::unmoving + Velocity<World3>({-8 * Metre / Second,
                                                   52 * Metre / Second,
                                                   0 * Metre / Second}), 0)));

  auto const conformal_map = similar_motion.conformal_map();

  auto const angular_velocity_of_from_frame =
      similar_motion.angular_velocity_of<World1>();
  auto const angular_velocity_of_to_frame =
      similar_motion.angular_velocity_of<World3>();

  auto const velocity_of_origin_of_from_frame =
      similar_motion.velocity_of_origin_of<World1>();
  auto const velocity_of_origin_of_to_frame =
      similar_motion.velocity_of_origin_of<World3>();

  auto const inverse = similar_motion.Inverse();

  auto const identity = SimilarMotion<World1, World3>::Identity();

  auto const q = World1::origin +
                 Displacement<World1>({-5 * Metre, 7 * Metre, 11 * Metre});
  auto const v = Velocity<World1>(
      {-2 * Metre / Second, 4 * Metre / Second, -8 * Metre / Second});

  EXPECT_THAT((similar_motion.Inverse() * similar_motion)({q, v}),
              Componentwise(AlmostEquals(q, 0), AlmostEquals(v, 0)));
  EXPECT_THAT(similar_motion.Inverse()(similar_motion({q, v})),
              Componentwise(AlmostEquals(q, 0), AlmostEquals(v, 0)));
}

}  // namespace physics
}  // namespace principia
