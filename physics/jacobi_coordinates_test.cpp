#include "jacobi_coordinates.hpp"

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using quantities::si::Kilogram;
using testing_utilities::AlmostEquals;
using testing_utilities::VanishesBefore;
using ::testing::ElementsAre;

namespace physics {

class JacobiCoordinatesTest : public ::testing::Test {
 protected:
  using Frame = geometry::Frame<serialization::Frame::TestTag,
                                serialization::Frame::TEST,
                                /*frame_is_inertial=*/true>;

  MassiveBody m1_ = MassiveBody(1 * Kilogram);
  MassiveBody m2_ = MassiveBody(2 * Kilogram);
};

TEST_F(JacobiCoordinatesTest, Jacobi) {
  auto const x_positions = [](JacobiCoordinates<Frame> const& system) {
    std::vector<Length> result;
    auto const barycentric_coordinates = system.BarycentricCoordinates();
    std::transform(barycentric_coordinates.begin(),
                   barycentric_coordinates.end(),
                   std::back_inserter(result),
                   [](RelativeDegreesOfFreedom<Frame> const& dof) {
                     return dof.displacement().coordinates().x;
                   });
    return result;
  };

  // e, i, ?, ?, and mean anomaly are 0.
  KeplerianElements<Frame> elements;

  JacobiCoordinates<Frame> system(m2_);
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
      ElementsAre(-1.0 / 3.0 * Metre, AlmostEquals(2.0 / 3.0 * Metre, 1)));

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
              ElementsAre(AlmostEquals(-2 * Metre, 1),
                          AlmostEquals(-1 * Metre, 2),
                          VanishesBefore(1 * Metre, 1),
                          5 * Metre));
}

}  // namespace physics
}  // namespace principia
