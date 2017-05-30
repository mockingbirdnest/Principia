
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
  Rotation<World, Camera> const world_to_camera_rotation(
      π / 6 * Radian,
      π / 4 * Radian,
      π / 3 * Radian,
      CardanoAngles::ZYX,
      DefinesFrame<Camera>());
  AffineMap<World, Camera, Length, OrthogonalMap> const world_to_camera_affine(
      World::origin + Displacement<World>({1 * Metre, 2 * Metre, -3 * Metre}),
      Camera::origin,
      world_to_camera_rotation.Forget());
  Perspective<World, Camera, Length, OrthogonalMap> perspective(
      world_to_camera_affine,
      /*focal=*/10 * Metre);

  Point<Displacement<World>> const point =
      World::origin + Displacement<World>({5 * Metre, 7 * Metre, 11 * Metre});
  EXPECT_THAT(
      world_to_camera_rotation(
          Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre})),
      Componentwise(
          AlmostEquals(Sqrt(3.0) / Sqrt(8.0) * Metre, 1),
          AlmostEquals((3.0 / Sqrt(32.0) - 1.0 / 4.0) * Metre, 0),
          AlmostEquals(Sqrt(3.0) / 4.0 * (1.0 + 1.0 / Sqrt(2.0)) * Metre, 1)));
  EXPECT_THAT(
      world_to_camera_rotation(
          Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre})),
      Componentwise(
          AlmostEquals(1.0 / Sqrt(8.0) * Metre, 0),
          AlmostEquals(Sqrt(3.0) / 4.0 * (1.0 + 1.0 / Sqrt(2.0)) * Metre, 2),
          AlmostEquals((1.0 / Sqrt(32.0) - 3.0 / 4.0) * Metre, 1)));
  EXPECT_THAT(world_to_camera_rotation(
                  Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre})),
              Componentwise(AlmostEquals(-1.0 / Sqrt(2.0) * Metre, 1),
                            AlmostEquals(Sqrt(3.0) / Sqrt(8.0) * Metre, 0),
                            AlmostEquals(1.0 / Sqrt(8.0) * Metre, 4)));
  EXPECT_THAT(perspective(point),
              Eq(RP2Element<Length>(1 * Metre, 1 * Metre, 1)));
}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
