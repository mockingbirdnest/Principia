
#include "jacobi_coordinates.hpp"

#include <map>
#include <vector>

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

  // e, i, Ω, ω, and mean anomaly are 0.
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

TEST_F(JacobiCoordinatesTest, Hierarchical) {

  // e, i, Ω, ω, and mean anomaly are 0.
  KeplerianElements<Frame> elements;

  // Invariant: |body_indices[bodies[i]] == i| for all |i|.
  std::map<not_null<MassiveBody const*>, int> body_indices;
  std::vector<not_null<MassiveBody const*>> bodies;

  auto const new_body = [&body_indices, &bodies](Mass const& mass) {
    auto body = make_not_null_unique<MassiveBody>(mass);
    bodies.emplace_back(body.get());
    body_indices[body.get()] = body_indices.size();
    return body;
  };

  // We construct a system as follows, where the |body_indices|
  // are, from left to right, 0, 2, 3, 1.
  // |<1 m>|     |<1 m>|
  // 2     1     1     2
  //   |<   7/3 m   >|
  std::vector<int> left_to_right = {0, 2, 3, 1};

  HierarchicalSystem<Frame> system(new_body(2 * Kilogram));
  elements.semimajor_axis = 7.0 / 3.0 * Metre;
  system.Add(new_body(2 * Kilogram), /*parent=*/bodies[0], elements);
  elements.semimajor_axis = 1 * Metre;
  system.Add(new_body(1 * Kilogram), /*parent=*/bodies[0], elements);
  elements.mean_anomaly = π * Radian;
  system.Add(new_body(1 * Kilogram), /*parent=*/bodies[1], elements);

  auto barycentric_system = system.Get();
  std::vector<int> expected_order = {0, 2, 1, 3};
  for (int i = 0; i < barycentric_system.bodies.size(); ++i) {
    LOG(ERROR) << body_indices[barycentric_system.bodies[i].get()];
    //EXPECT_TRUE(bodies[i] == barycentric_system.bodies[i].get()) << i;
  }
  std::vector<Length> x_positions;
  std::transform(barycentric_system.degrees_of_freedom.begin(),
                 barycentric_system.degrees_of_freedom.end(),
                 std::back_inserter(x_positions),
                 [](DegreesOfFreedom<Frame> const& dof) {
                   return (dof.position() - Frame::origin).coordinates().x;
                 });
  EXPECT_THAT(x_positions,
              ElementsAre(AlmostEquals(-1.5 * Metre, 1),
                          AlmostEquals(-0.5 * Metre, 2),
                          AlmostEquals(0.5 * Metre, 1),
                          AlmostEquals(1.5 * Metre, 2)));
}

}  // namespace physics
}  // namespace principia
