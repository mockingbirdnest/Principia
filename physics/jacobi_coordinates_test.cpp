#include "physics/jacobi_coordinates.hpp"

#include <algorithm>
#include <vector>

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using ::testing::ElementsAre;
using namespace principia::geometry::_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_jacobi_coordinates;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_vanishes_before;

class JacobiCoordinatesTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag, Inertial>;

  MassiveBody m1_ = MassiveBody(1 * Kilogram);
  MassiveBody m2_ = MassiveBody(2 * Kilogram);
};

TEST_F(JacobiCoordinatesTest, Jacobi) {
  auto const x_positions = [](JacobiCoordinates<World> const& system) {
    std::vector<Length> result;
    auto const barycentric_dof = system.BarycentricDegreesOfFreedom();
    std::transform(barycentric_dof.begin(),
                   barycentric_dof.end(),
                   std::back_inserter(result),
                   [](RelativeDegreesOfFreedom<World> const& dof) {
                     return dof.displacement().coordinates().x;
                   });
    return result;
  };

  // i, and Ω are 0 by default.
  KeplerianElements<World> elements;
  elements.eccentricity = 0;
  elements.argument_of_periapsis = 0 * Radian;
  elements.mean_anomaly = 0 * Radian;

  JacobiCoordinates<World> system(m2_);
  EXPECT_EQ(2 * Kilogram, system.System().mass());
  EXPECT_THAT(x_positions(system), ElementsAre(0 * Metre));

  elements.semimajor_axis = 1 * Metre;
  system.Add(m1_, elements);
  // The system now consists of a 2 kg mass and a 1 kg mass, with the barycentre
  // one third of the way, as shown.
  // 2  1
  //  ^ barycentre
  EXPECT_EQ(3 * Kilogram, system.System().mass());
  EXPECT_THAT(
      x_positions(system),
      ElementsAre(AlmostEquals(-1.0 / 3.0 * Metre, 0),
                  AlmostEquals(2.0 / 3.0 * Metre, 1)));

  elements.semimajor_axis = 5.0 / 3.0 * Metre;
  system.Add(m2_, elements);
  // 2  1  2
  //    ^ barycentre
  EXPECT_EQ(5 * Kilogram, system.System().mass());
  EXPECT_THAT(x_positions(system),
              ElementsAre(-1 * Metre, 0 * Metre, 1 * Metre));

  elements.semimajor_axis = 6 * Metre;
  system.Add(m1_, elements);
  // 2  1  2  .  .  .  .  1
  //       ^ barycentre
  EXPECT_THAT(x_positions(system),
              ElementsAre(AlmostEquals(-2 * Metre, 0),
                          AlmostEquals(-1 * Metre, 0),
                          VanishesBefore(1 * Metre, 0),
                          5 * Metre));
}

}  // namespace physics
}  // namespace principia
