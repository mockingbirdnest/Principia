
#include "jacobi_coordinates.hpp"

#include <algorithm>
#include <vector>

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_jacobi_coordinates {

using geometry::Frame;
using geometry::Inertial;
using quantities::Length;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using testing_utilities::AlmostEquals;
using testing_utilities::VanishesBefore;
using ::testing::ElementsAre;

class JacobiCoordinatesTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag, Inertial>;

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
      ElementsAre(AlmostEquals(-1.0 / 3.0 * Metre, 1),
                  AlmostEquals(2.0 / 3.0 * Metre, 0)));

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

}  // namespace internal_jacobi_coordinates
}  // namespace physics
}  // namespace principia
