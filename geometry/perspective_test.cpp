
#include <limits>

#include "geometry/affine_map.hpp"
#include "geometry/frame.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "geometry/rp2_point.hpp"
#include "geometry/sphere.hpp"
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
using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::SizeIs;
using ::testing::_;

class PerspectiveTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST1>;
  using Camera = Frame<serialization::Frame::TestTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::TEST2>;
};

TEST_F(PerspectiveTest, Basic) {
  Position<World> const camera_origin =
      World::origin + Displacement<World>({1 * Metre, 2 * Metre, -3 * Metre});
  Rotation<World, Camera> const world_to_camera_rotation(
      π / 6 * Radian,
      π / 4 * Radian,
      π / 3 * Radian,
      CardanoAngles::ZYX,
      DefinesFrame<Camera>());
  RigidTransformation<World, Camera> const world_to_camera_transformation(
      camera_origin,
      Camera::origin,
      world_to_camera_rotation.Forget<OrthogonalMap>());
  Perspective<World, Camera> perspective(world_to_camera_transformation,
                                         /*focal=*/10 * Metre);

  // Check that points in the camera z axis get projected to the origin of ℝP².
  Displacement<World> const camera_z_axis = world_to_camera_rotation.Inverse()(
      Displacement<Camera>({0 * Metre, 0 * Metre, 1 * Metre}));
  Position<World> const p0 = camera_origin;
  Position<World> const p1 = camera_origin + 1 * camera_z_axis;
  Position<World> const p2 = camera_origin + 10 * camera_z_axis;
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
  Position<World> const p3 = p1 + 5 * camera_x_axis;
  Position<World> const p4 = p1 + 7 * camera_x_axis;
  EXPECT_THAT(perspective(p3).y(), VanishesBefore(1 * Metre, 20));
  EXPECT_THAT(perspective(p4).y(), VanishesBefore(1 * Metre, 10));

  // Check that points on the camera y axis get projected on the y axis of ℝP².
  Displacement<World> const camera_y_axis = world_to_camera_rotation.Inverse()(
      Displacement<Camera>({0 * Metre, 1 * Metre, 0 * Metre}));
  Position<World> const p5 = p1 - 11 * camera_y_axis;
  Position<World> const p6 = p1 + 13 * camera_y_axis;
  EXPECT_THAT(perspective(p5).x(), VanishesBefore(1 * Metre, 120));
  EXPECT_THAT(perspective(p6).x(), VanishesBefore(1 * Metre, 0));

  // Check that aligned points are aligned in ℝP².
  Position<World> const p7 =
      camera_origin +
      Displacement<World>({17 * Metre, -23 * Metre, 29 * Metre});
  Position<World> const p8 =
      camera_origin +
      Displacement<World>({18 * Metre, -21 * Metre, 24 * Metre});
  Position<World> const p9 =
      camera_origin +
      Displacement<World>({19 * Metre, -19 * Metre, 19 * Metre});
  auto const q7 = perspective(p7);
  auto const q8 = perspective(p8);
  auto const q9 = perspective(p9);
  EXPECT_THAT((q8.x() - q7.x()) * (q9.y() - q7.y()) -
                  (q9.x() - q7.x()) * (q8.y() - q7.y()),
              VanishesBefore(1 * Metre * Metre, 6));

  // Check that the focal works as expected.
  Position<World> const p10 =
      camera_origin + 1 * camera_x_axis + 2 * camera_y_axis + 3 * camera_z_axis;
  EXPECT_THAT(perspective(p10),
              Componentwise(AlmostEquals(1.0 / 0.3 * Metre, 3),
                            AlmostEquals(2.0 / 0.3 * Metre, 2)));
}

TEST_F(PerspectiveTest, SegmentBehindFocalPlane) {
  Perspective<World, Camera> perspective(
      RigidTransformation<World, Camera>::Identity(),
      /*focal=*/1 * Metre);

  // In front of the camera.
  Position<World> const p1 =
      World::origin +
      Displacement<World>({1 * Metre, 2 * Metre, 3 * Metre});
  // Behind the camera.
  Position<World> const p2 =
      World::origin +
      Displacement<World>({4 * Metre, 5 * Metre, -6 * Metre});

  auto const segment = perspective.SegmentBehindFocalPlane({p1, p2});
  EXPECT_TRUE(segment.has_value());
  EXPECT_THAT(segment->first, AlmostEquals(p1, 0));
  EXPECT_THAT(
      segment->second,
      AlmostEquals(World::origin +
                   Displacement<World>(
                       {5.0 / 3.0 * Metre, 8.0 / 3.0 * Metre, 1 * Metre}),
                   1));
}

TEST_F(PerspectiveTest, IsHiddenBySphere) {
  Perspective<World, Camera> perspective(
      RigidTransformation<World, Camera>::Identity(),
      /*focal=*/1 * Metre);

  Sphere<World> const sphere(
      World::origin + Displacement<World>({10 * Metre, 20 * Metre, 30 * Metre}),
      /*radius=*/3 * Metre);

  // Within the sphere.
  Position<World> const p1 =
      World::origin +
      Displacement<World>({11 * Metre, 19 * Metre, 32 * Metre});
  // Far from the sphere.
  Position<World> const p2 =
      World::origin +
      Displacement<World>({100 * Metre, 50 * Metre, -70 * Metre});
  // Behind the sphere.
  Position<World> const p3 =
      World::origin +
      Displacement<World>({100 * Metre, 202 * Metre, 305 * Metre});
  // In front of the sphere.
  Position<World> const p4 =
      World::origin +
      Displacement<World>({2 * Metre, 4.05 * Metre, 6 * Metre});

  EXPECT_TRUE(perspective.IsHiddenBySphere(p1, sphere));
  EXPECT_FALSE(perspective.IsHiddenBySphere(p2, sphere));
  EXPECT_TRUE(perspective.IsHiddenBySphere(p3, sphere));
  EXPECT_FALSE(perspective.IsHiddenBySphere(p4, sphere));
}

TEST_F(PerspectiveTest, SphereSin²HalfAngle) {
  Perspective<World, Camera> perspective(
      RigidTransformation<World, Camera>::Identity(),
      /*focal=*/1 * Metre);

  Sphere<World> const sphere(
      World::origin + Displacement<World>({0 * Metre, 0 * Metre, 100 * Metre}),
      /*radius=*/1 * Metre);

  EXPECT_THAT(perspective.SphereSin²HalfAngle(sphere),
              AlmostEquals(0.0001, 0));
}

TEST_F(PerspectiveTest, Output) {
  Perspective<World, Camera> perspective(
      RigidTransformation<World, Camera>::Identity(),
      /*focal=*/1 * Metre);
  std::cout << perspective << "\n";
}

class VisibleSegmentsTest : public PerspectiveTest {
 protected:
  VisibleSegmentsTest()
      :  // The camera is on the x-axis and looks towards the positive x.
        camera_origin_(
            World::origin +
            Displacement<World>({-10 * Metre, 0 * Metre, 0 * Metre})),
        world_to_camera_transformation_(
            camera_origin_,
            Camera::origin,
            OrthogonalMap<World, Camera>::Identity()),
        perspective_(world_to_camera_transformation_,
                     /*focal=*/1 * Metre),
        // The sphere is at the origin and has unit radius.
        sphere_(World::origin,
                /*radius=*/1 * Metre) {
    // The equation of the line KP in the plane x-z is:
    //   x - 3 * sqrt(11) z + 10 = 0
    // This is used to compute many of the expected values below.
  }

  Position<World> const camera_origin_;
  RigidTransformation<World, Camera> const world_to_camera_transformation_;
  Perspective<World, Camera> const perspective_;
  Sphere<World> const sphere_;
};

// A segment away from the sphere with x > 0.
TEST_F(VisibleSegmentsTest, AwayPositiveX) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({10 * Metre, 20 * Metre, 30 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({9 * Metre, 21 * Metre, 32 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// A segment away from the sphere with x < 0.
TEST_F(VisibleSegmentsTest, AwayNegativeX) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-5 * Metre, 20 * Metre, 30 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({-3 * Metre, 21 * Metre, 32 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// A segment tangent to the sphere when seen from the camera.  Need a 1-ulp
// tolerance to make sure that the segment doesn't get cut in two because of
// rounding.
TEST_F(VisibleSegmentsTest, TangentNotBitten) {
  double const ε = std::numeric_limits<double>::epsilon();
  Position<World> const p1 =
      World::origin +
      Displacement<World>(
          {9.8 * Metre, Sqrt(3.96) * (1 + ε) * Metre, 7 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>(
          {9.8 * Metre, Sqrt(3.96) * (1 + ε) * Metre, -9 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// Same as above but the segment gets cut in two.
TEST_F(VisibleSegmentsTest, TangentBittenBittenBitten) {
  double const ε = std::numeric_limits<double>::epsilon();
  Position<World> const p1 =
      World::origin +
      Displacement<World>(
          {9.8 * Metre, Sqrt(3.96) * (1 - ε) * Metre, 7 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>(
          {9.8 * Metre, Sqrt(3.96) * (1 - ε) * Metre, -9 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(Pair(p1, _), Pair(_, p2)));
}

// A segment entirely in front of the sphere, smaller than the sphere.
TEST_F(VisibleSegmentsTest, InFrontSmaller) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-5 * Metre, 0 * Metre, 0.1 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({-5 * Metre, 0 * Metre, -0.1 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// A segment entirely in front of the sphere, out of the sphere in one
// direction.
TEST_F(VisibleSegmentsTest, InFrontOnTheSide) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-5 * Metre, 0 * Metre, 3 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({-5 * Metre, 0 * Metre, -0.1 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// A segment entirely in front of the sphere, out of the sphere in both
// directions.
TEST_F(VisibleSegmentsTest, InFrontLarger) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-5 * Metre, 0 * Metre, 3 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({-5 * Metre, 0 * Metre, -5 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// A segment entirely in front of the sphere, not parallel to the z-axis.
TEST_F(VisibleSegmentsTest, InFrontNotAlongZ) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-5 * Metre, 0.2 * Metre, -3 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({-5 * Metre, 3 * Metre, 0.1 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// A segment behind the sphere, with a length smaller than the diametre.
TEST_F(VisibleSegmentsTest, BehindSmallerThanDiameter) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({5 * Metre, 0.9 * Metre, 0.1 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({5 * Metre, -0.8 * Metre, -0.3 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_), IsEmpty());
}

// A segment behind the sphere, with a length larger than the diametre.
TEST_F(VisibleSegmentsTest, BehindLargerThanDiameter) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({10 * Metre, 1.5 * Metre, -1.1 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({10 * Metre, -1.6 * Metre, 1.2 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_), IsEmpty());
}

// A segment intersecting the front of the sphere and not intersecting the
// cone.
TEST_F(VisibleSegmentsTest, IntersectingFrontOfTheSphere) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-0.5 * Metre, 0 * Metre, -1 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({-0.5 * Metre, 0 * Metre, 2 * Metre});
  Position<World> const p3 =
      World::origin +
      Displacement<World>({-0.5 * Metre, 0 * Metre, -Sqrt(3.0) / 2 * Metre});
  Position<World> const p4 =
      World::origin +
      Displacement<World>({-0.5 * Metre, 0 * Metre, Sqrt(3.0) / 2 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(Pair(p1, AlmostEquals(p3, 0)),
                          Pair(AlmostEquals(p4, 2), p2)));
}

// A segment intersecting the cone in front of the centre of the sphere, both
// extremities are visible.
TEST_F(VisibleSegmentsTest, IntersectingConeInFrontOfTheSphereCentre) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-0.05 * Metre, 0 * Metre, -3 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({-0.05 * Metre, 0 * Metre, 2 * Metre});
  Position<World> const p3 =
      World::origin +
      Displacement<World>(
          {-0.05 * Metre, 0 * Metre, -199.0 / (60.0 * Sqrt(11.0)) * Metre});
  Position<World> const p4 =
      World::origin +
      Displacement<World>(
          {-0.05 * Metre, 0 * Metre, 199.0 / (60.0 * Sqrt(11.0)) * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(Pair(p1, AlmostEquals(p3, 4)),
                          Pair(AlmostEquals(p4, 4), p2)));
}

// A segment intersecting the cone behind the centre of the sphere, both
// extremities are visible.
TEST_F(VisibleSegmentsTest, IntersectingConeTwoVisibleSegments) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({10 * Metre, 0 * Metre, -3 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({10 * Metre, 0 * Metre, 5 * Metre});
  Position<World> const p3 =
      World::origin +
      Displacement<World>(
          {10 * Metre, 0 * Metre, -20.0 / (3.0 * Sqrt(11.0)) * Metre});
  Position<World> const p4 =
      World::origin +
      Displacement<World>(
          {10 * Metre, 0 * Metre, 20.0 / (3.0 * Sqrt(11.0)) * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(Pair(p1, AlmostEquals(p3, 3)),
                          Pair(AlmostEquals(p4, 15), p2)));
}

// A segment intersecting the cone behind the centre of the sphere, only one
// extremity is visible.
TEST_F(VisibleSegmentsTest, IntersectingConeOneVisibleSegment) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({9 * Metre, 0 * Metre, 3 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>(
          {11 * Metre, 0 * Metre, (40.0 / (3.0 * Sqrt(11.0)) - 3.0) * Metre});
  Position<World> const p3 =
      World::origin +
      Displacement<World>(
          {10 * Metre, 0 * Metre, 20.0 / (3.0 * Sqrt(11.0)) * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(Pair(p1, AlmostEquals(p3, 95))));
}

// A segment intersecting the sphere on one side and the cone behind the
// centre of the sphere on the other side.  Both extremities are visible.
TEST_F(VisibleSegmentsTest, IntersectingConeAndSphere) {
  Position<World> const p1 =
      World::origin +
      Displacement<World>({-2 * Metre, 0 * Metre, -2 * Metre});
  Position<World> const p2 =
      World::origin +
      Displacement<World>({3 * Metre, 0 * Metre, 3 * Metre});
  Position<World> const p3 =
      World::origin +
      Displacement<World>(
          {-Sqrt(0.5) * Metre, 0 * Metre, -Sqrt(0.5) * Metre});
  Position<World> const p4 =
      World::origin +
      Displacement<World>({10.0 / (3.0 * Sqrt(11.0) - 1.0) * Metre,
                           0 * Metre,
                           10.0 / (3.0 * Sqrt(11.0) - 1.0) * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(Pair(p1, AlmostEquals(p3, 1)),
                          Pair(AlmostEquals(p4, 4), p2)));
}

// A segment vaguely parallel to the axis of the cone.  It does intersect the
// cone twice, but one of the intersections is behind the camera so there is no
// hiding.
TEST_F(VisibleSegmentsTest, HyperbolicIntersection) {
  Position<World> const p1 =
      World::origin + Displacement<World>({+7.77429986470929890e+00 * Metre,
                                           -6.96329812586472841e+00 * Metre,
                                           +4.69986567261152999e+00 * Metre});
  Position<World> const p2 =
      World::origin + Displacement<World>({-6.71731874077453561e+00 * Metre,
                                           -7.93922086105482183e+00 * Metre,
                                           +4.06349175030368848e+00 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, sphere_),
              ElementsAre(segment));
}

// A case where the intersections are far outside of the segment.  This used to
// be mishandled.
TEST_F(VisibleSegmentsTest, BehindCamera) {
  Position<World> const camera_origin(
      World::origin + Displacement<World>({-1.35994803226833153e+10 * Metre,
                                           +6.48711944107992947e+06 * Metre,
                                           +1.16940398868560791e+05 * Metre}));
  RigidTransformation<World, Camera> world_to_camera_transformation(
      camera_origin, Camera::origin, OrthogonalMap<World, Camera>::Identity());
  Perspective<World, Camera> const perspective(
      world_to_camera_transformation,
      /*focal=*/1.00000001794823179e+00 * Metre);
  Sphere<World> const sphere(
      World::origin + Displacement<World>({-1.35998847769040604e+10 * Metre,
                                           +5.65312885169138480e+06 * Metre,
                                           -2.45548475776178384e+03 * Metre}),
      /*radius=*/+6.30000000000000000e+05 * Metre);
  Position<World> const p1 =
      World::origin + Displacement<World>({-1.35993763102505913e+10 * Metre,
                                           +1.02319916576216407e+07 * Metre,
                                           -1.54814310483471490e+03 * Metre});
  Position<World> const p2 =
      World::origin + Displacement<World>({-1.35993763068824501e+10 * Metre,
                                           +1.02318443313949741e+07 * Metre,
                                           -1.54719462661686703e+03 * Metre});

  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective.VisibleSegments(segment, sphere),
              ElementsAre(segment));
}

// A case where the camera is inside the sphere.
TEST_F(VisibleSegmentsTest, CameraInsideSphere) {
  Position<World> const camera_origin(
      World::origin + Displacement<World>({-1.35993502776454182e+10 * Metre,
                                           +4.79792595486199576e+06 * Metre,
                                           -1.57980008827209473e+05 * Metre}));
  RigidTransformation<World, Camera> world_to_camera_transformation(
      camera_origin, Camera::origin, OrthogonalMap<World, Camera>::Identity());
  Perspective<World, Camera> const perspective(world_to_camera_transformation,
                                               /*focal=*/1 * Metre);
  Sphere<World> const sphere(
      World::origin + Displacement<World>({-1.35998843607365036e+10 * Metre,
                                           +5.00094497194097005e+06 * Metre,
                                           -2.45569064268199872e+03 * Metre}),
      /*radius=*/+6.30000000000000000e+05 * Metre);
  Position<World> const p1 =
      World::origin + Displacement<World>({-1.35993763102505913e+10 * Metre,
                                           +1.02319916576216407e+07 * Metre,
                                           -1.54814310483471490e+03 * Metre});
  Position<World> const p2 =
      World::origin + Displacement<World>({-1.35993763068824501e+10 * Metre,
                                           +1.02318443313949741e+07 * Metre,
                                           -1.54719462661686703e+03 * Metre});

  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective.VisibleSegments(segment, sphere), IsEmpty());
}

// A hyperbolic case where the segment is entirely hidden.  It used to be
// entirely visible.
TEST_F(VisibleSegmentsTest, AnotherHyperbolicIntersection) {
  Position<World> const camera_origin(
      World::origin + Displacement<World>({-1.35999109873531647e+10 * Metre,
                                           -2.00764947838563850e+05 * Metre,
                                           -2.22307361124038696e+05 * Metre}));
  RigidTransformation<World, Camera> world_to_camera_transformation(
      camera_origin, Camera::origin, OrthogonalMap<World, Camera>::Identity());
  Perspective<World, Camera> const perspective(world_to_camera_transformation,
                                               /*focal=*/1 * Metre);
  Sphere<World> const sphere(
      World::origin + Displacement<World>({-1.35998825622369633e+10 * Metre,
                                           +2.59904061185893603e+06 * Metre,
                                           -2.45644535030410043e+03 * Metre}),
      /*radius=*/+6.30000000000000000e+05 * Metre);
  Position<World> const p1 =
      World::origin + Displacement<World>({-1.35993763102505913e+10 * Metre,
                                           +1.02319916576216407e+07 * Metre,
                                           -1.54814310483471490e+03 * Metre});
  Position<World> const p2 =
      World::origin + Displacement<World>({-1.35993763068824501e+10 * Metre,
                                           +1.02318443313949741e+07 * Metre,
                                           -1.54719462661686703e+03 * Metre});

  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective.VisibleSegments(segment, sphere), IsEmpty());
}

// A case where the segment is seen under a very small angle (3.2e-08 radian)
// from the camera.
TEST_F(VisibleSegmentsTest, SmallAngle) {
  Position<World> const camera_origin =
      World::origin + Displacement<World>({+6.91404109651565552e+05 * Metre,
                                           -6.24300985912289470e+04 * Metre,
                                           -1.05537869332772680e+04 * Metre});
  RigidTransformation<World, Camera> world_to_camera_transformation(
      camera_origin, Camera::origin, OrthogonalMap<World, Camera>::Identity());
  Perspective<World, Camera> const perspective(world_to_camera_transformation,
                                               /*focal=*/1 * Metre);
  Sphere<World> const sphere(
      World::origin + Displacement<World>({+0.00000000000000000e+00 * Metre,
                                           +0.00000000000000000e+00 * Metre,
                                           +0.00000000000000000e+00 * Metre}),
      /*radius=*/+6.30000000000000000e+05 * Metre);
  Position<World> const p1 =
      World::origin + Displacement<World>({+6.66879854196548462e+05 * Metre,
                                           +7.60680419044313021e+05 * Metre,
                                           +1.90194540126288375e+04 * Metre});
  Position<World> const p2 =
      World::origin + Displacement<World>({+6.66878618906021118e+05 * Metre,
                                           +7.60722757674634922e+05 * Metre,
                                           +1.90209760584930445e+04 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective.VisibleSegments(segment, sphere),
              ElementsAre(segment));
}

TEST_F(VisibleSegmentsTest, MultipleSpheres) {
  Sphere<World> const sphere2(
    World::origin + Displacement<World>({0 * Metre, 0 * Metre, 5 * Metre}),
      /*radius=*/1 * Metre);
  Position<World> const p1 =
      World::origin + Displacement<World>({2 * Metre, 0 * Metre, -10 * Metre});
  Position<World> const p2 =
      World::origin + Displacement<World>({2 * Metre, 0 * Metre, 10 * Metre});
  Segment<World> segment{p1, p2};
  EXPECT_THAT(perspective_.VisibleSegments(segment, {sphere_, sphere2}),
              SizeIs(3));
}

}  // namespace internal_perspective
}  // namespace geometry
}  // namespace principia
