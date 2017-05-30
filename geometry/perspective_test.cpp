
#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/frame.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "geometry/rp2_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace geometry {
namespace internal_perspective {

using quantities::Length;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Radian;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
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
  EXPECT_EQ(RP2Element<Length>(0 * Metre, 0 * Metre, 1), perspective(p1));
  EXPECT_EQ(RP2Element<Length>(0 * Metre, 0 * Metre, 1), perspective(p2));

  // Check that points on the camera x axis get projected on the x axis of ℝP².
  Displacement<World> const camera_x_axis = world_to_camera_rotation.Inverse()(
      Displacement<Camera>({1 * Metre, 0 * Metre, 0 * Metre}));
  Point<Displacement<World>> const p3 = camera_origin + 5 * camera_x_axis;
  Point<Displacement<World>> const p4 = camera_origin + 7 * camera_x_axis;
  EXPECT_EQ(0 * Metre, perspective(p3).y());
  EXPECT_EQ(0 * Metre, perspective(p4).y());

  // Check that points on the camera y axis get projected on the y axis of ℝP².
  Displacement<World> const camera_y_axis = world_to_camera_rotation.Inverse()(
      Displacement<Camera>({0 * Metre, 1 * Metre, 0 * Metre}));
  Point<Displacement<World>> const p5 = camera_origin - 11 * camera_y_axis;
  Point<Displacement<World>> const p6 = camera_origin + 13 * camera_y_axis;
  EXPECT_EQ(0 * Metre, perspective(p5).x());
  EXPECT_EQ(0 * Metre, perspective(p6).x());

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
  RP2Element<Length> const q7 = perspective(p7);
  RP2Element<Length> const q8 = perspective(p8);
  RP2Element<Length> const q9 = perspective(p9);
  EXPECT_EQ(0 * Metre * Metre,
            (q8.x() - q7.x()) * (q9.y() - q7.y()) -
                (q9.x() - q7.x()) * (q8.y() - q7.y()));
}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
