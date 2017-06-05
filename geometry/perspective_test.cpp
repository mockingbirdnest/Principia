
#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/frame.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "geometry/rp2_point.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

using quantities::Length;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Radian;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;
using ::testing::Eq;

class PerspectiveTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, false>;
  using Camera = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST2, false>;
};

TEST_F(PerspectiveTest, Basic) {
  Point<Displacement<World>> const camera_origin =
      World::origin + Displacement<World>({1 * Metre, 2 * Metre, -3 * Metre});
  Rotation<World, Camera> const world_to_camera_rotation(
      π / 6 * Radian,
      π / 4 * Radian,
      π / 3 * Radian,
      CardanoAngles::ZYX,
      DefinesFrame<Camera>());
  AffineMap<World, Camera, Length, OrthogonalMap> const world_to_camera_affine(
      camera_origin,
      Camera::origin,
      world_to_camera_rotation.Forget());
  Perspective<World, Camera, Length, OrthogonalMap> perspective(
      world_to_camera_affine,
      /*focal=*/10 * Metre);

  // Check that points in the camera z axis get projected to the origin of ℝP².
  Displacement<World> const camera_z_axis = world_to_camera_rotation.Inverse()(
      Displacement<Camera>({0 * Metre, 0 * Metre, 1 * Metre}));
  Point<Displacement<World>> const p0 = camera_origin;
  Point<Displacement<World>> const p1 = camera_origin + 1 * camera_z_axis;
  Point<Displacement<World>> const p2 = camera_origin + 10 * camera_z_axis;
  EXPECT_TRUE(perspective(p0).is_at_infinity());
  EXPECT_THAT(perspective(p1),
              Componentwise(VanishesBefore(1 * Metre, 10),
                            VanishesBefore(1 * Metre, 0)));
  EXPECT_THAT(perspective(p2),
              Componentwise(VanishesBefore(1 * Metre, 8),
                            VanishesBefore(1 * Metre, 4)));

  // Check that points on the camera x axis get projected on the x axis of ℝP².
  Displacement<World> const camera_x_axis = world_to_camera_rotation.Inverse()(
      Displacement<Camera>({1 * Metre, 0 * Metre, 0 * Metre}));
  Point<Displacement<World>> const p3 = p1 + 5 * camera_x_axis;
  Point<Displacement<World>> const p4 = p1 + 7 * camera_x_axis;
  EXPECT_THAT(perspective(p3).y(), VanishesBefore(1 * Metre, 20));
  EXPECT_THAT(perspective(p4).y(), VanishesBefore(1 * Metre, 10));

  // Check that points on the camera y axis get projected on the y axis of ℝP².
  Displacement<World> const camera_y_axis = world_to_camera_rotation.Inverse()(
      Displacement<Camera>({0 * Metre, 1 * Metre, 0 * Metre}));
  Point<Displacement<World>> const p5 = p1 - 11 * camera_y_axis;
  Point<Displacement<World>> const p6 = p1 + 13 * camera_y_axis;
  EXPECT_THAT(perspective(p5).x(), VanishesBefore(1 * Metre, 120));
  EXPECT_THAT(perspective(p6).x(), VanishesBefore(1 * Metre, 0));

  // Check that aligned points are aligned in ℝP².
  Point<Displacement<World>> const p7 =
      camera_origin +
      Displacement<World>({17 * Metre, -23 * Metre, 29 * Metre});
  Point<Displacement<World>> const p8 =
      camera_origin +
      Displacement<World>({18 * Metre, -21 * Metre, 24 * Metre});
  Point<Displacement<World>> const p9 =
      camera_origin +
      Displacement<World>({19 * Metre, -19 * Metre, 19 * Metre});
  RP2Point<Length, Camera> const q7 = perspective(p7);
  RP2Point<Length, Camera> const q8 = perspective(p8);
  RP2Point<Length, Camera> const q9 = perspective(p9);
  EXPECT_THAT((q8.x() - q7.x()) * (q9.y() - q7.y()) -
                  (q9.x() - q7.x()) * (q8.y() - q7.y()),
              VanishesBefore(1 * Metre * Metre, 6));

  // Check that the focal works as expected.
  Point<Displacement<World>> const p10 =
      camera_origin + 1 * camera_x_axis + 2 * camera_y_axis + 3 * camera_z_axis;
  EXPECT_THAT(perspective(p10),
              Componentwise(AlmostEquals(1.0 / 0.3 * Metre, 3),
                            AlmostEquals(2.0 / 0.3 * Metre, 2)));
}

TEST_F(PerspectiveTest, IsHiddenBySphere) {
  Perspective<World, Camera, Length, OrthogonalMap> perspective(
      AffineMap<World, Camera, Length, OrthogonalMap>::Identity(),
      /*focal=*/1 * Metre);

  Point<Displacement<World>> const centre =
      World::origin +
      Displacement<World>({10 * Metre, 20 * Metre, 30 * Metre});
  Length const radius = 3 * Metre;

  // Within the sphere.
  Point<Displacement<World>> const p1 =
      centre +
      Displacement<World>({1 * Metre, -1 * Metre, 2 * Metre});
  // Far from the sphere.
  Point<Displacement<World>> const p2 =
      World::origin +
      Displacement<World>({100 * Metre, 50 * Metre, -70 * Metre});
  // Behind the sphere.
  Point<Displacement<World>> const p3 =
      World::origin +
      Displacement<World>({100 * Metre, 202 * Metre, 305 * Metre});
  // In front of the sphere.
  Point<Displacement<World>> const p4 =
      World::origin +
      Displacement<World>({2 * Metre, 5 * Metre, 6 * Metre});

  EXPECT_TRUE(perspective.IsHiddenBySphere(p1, centre, radius));
  EXPECT_FALSE(perspective.IsHiddenBySphere(p2, centre, radius));
  EXPECT_TRUE(perspective.IsHiddenBySphere(p3, centre, radius));
  EXPECT_FALSE(perspective.IsHiddenBySphere(p4, centre, radius));
}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
