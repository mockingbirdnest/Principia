
#include "ksp_plugin/planetarium.hpp"

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "gtest/gtest.h"
#include "physics/mock_dynamic_frame.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/numbers.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Displacement;
using geometry::LinearMap;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using physics::MockDynamicFrame;
using physics::RigidMotion;
using physics::RigidTransformation;
using quantities::Cos;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::VanishesBefore;
using ::testing::_;
using ::testing::Return;

class PlanetariumTest : public ::testing::Test {};

TEST_F(PlanetariumTest, PlotMethod0) {
  Sphere<Length, Barycentric> sphere(Barycentric::origin, /*radius=*/1 * Metre);

  // A circular trajectory around the origin.
  DiscreteTrajectory<Barycentric> trajectory;
  for (Time t; t <= 10 * Second; t += 1 * Second) {
    DegreesOfFreedom<Barycentric> const degrees_of_freedom(
        Barycentric::origin +
            Displacement<Barycentric>(
                {10 * Metre * Sin(2 * π * t * Radian / (10 * Second)),
                 10 * Metre * Cos(2 * π * t * Radian / (10 * Second)),
                 0 * Metre}),
        Velocity<Barycentric>());
    trajectory.Append(Instant() + t, degrees_of_freedom);
  }

  // The camera is located as {0, 20, 0} and is looking along -y.
  Perspective<Navigation, Camera, Length, OrthogonalMap> const perspective(
      AffineMap<Navigation, Camera, Length, OrthogonalMap>(
          Navigation::origin +
              Displacement<Navigation>({0 * Metre, 20 * Metre, 0 * Metre}),
          Camera::origin,
          Rotation<Navigation, Camera>(
              Vector<double, Navigation>({1, 0, 0}),
              Vector<double, Navigation>({0, 0, 1}),
              Bivector<double, Navigation>({0, -1, 0})).Forget()),
      /*focal=*/5 * Metre);

  MockDynamicFrame<Barycentric, Navigation> plotting_frame;
  EXPECT_CALL(plotting_frame, ToThisFrameAtTime(_))
      .WillRepeatedly(Return(RigidMotion<Barycentric, Navigation>(
          RigidTransformation<Barycentric, Navigation>::Identity(),
          AngularVelocity<Barycentric>(),
          Velocity<Barycentric>())));

  Planetarium planetarium({sphere}, perspective, &plotting_frame);
  auto const rp2_points =
      planetarium.PlotMethod0(trajectory, Instant() + 10 * Second);

  for (auto const& rp2_point : rp2_points) {
    // The following limit is obtained by elementary geometry by noticing that
    // the circle is viewed from the camera under an angle of π / 6.
    EXPECT_LE(rp2_point.x(), 5.0 / Sqrt(3.0) * Metre);
    EXPECT_THAT(rp2_point.y(), VanishesBefore(1 * Metre, 6, 13));
  }
}

}  // namespace internal_planetarium
}  // namespace ksp_plugin
}  // namespace principia
