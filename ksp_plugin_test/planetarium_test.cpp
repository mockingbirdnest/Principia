
#include "ksp_plugin/planetarium.hpp"

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "gtest/gtest.h"
#include "physics/mock_dynamic_frame.hpp"
#include "quantities/numbers.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::AffineMap;
using geometry::Bivector;
using geometry::Displacement;
using geometry::LinearMap;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using physics::MockDynamicFrame;
using quantities::Cos;
using quantities::Sin;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

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
  // TODO(phl): Why Bivector below?
  Perspective<Barycentric, Camera, Length, OrthogonalMap> const perspective(
      AffineMap<Barycentric, Camera, Length, OrthogonalMap>(
          Barycentric::origin +
              Displacement<Barycentric>({0 * Metre, 20 * Metre, 0 * Metre}),
          Camera::origin,
          Rotation<Barycentric, Camera>(
              Vector<double, Barycentric>({1, 0, 0}),
              Vector<double, Barycentric>({0, 0, 1}),
              Bivector<double, Barycentric>({0, -1, 0})).Forget()),
      /*focal=*/5 * Metre);

  MockDynamicFrame<Barycentric, Navigation> plotting_frame;
  Planetarium planetarium({sphere}, perspective, &plotting_frame);

  auto const rp2_points =
      planetarium.PlotMethod0(trajectory, Instant() + 10 * Second);
}

}  // namespace internal_planetarium
}  // namespace ksp_plugin
}  // namespace principia
