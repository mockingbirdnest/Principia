
#include "hierarchical_system.hpp"

#include <algorithm>
#include <map>
#include <vector>

#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_hierarchical_system {

using base::make_not_null_unique;
using geometry::Frame;
using geometry::Inertial;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Sqrt;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::VanishesBefore;
using ::testing::ElementsAre;

class HierarchicalSystemTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag, Inertial>;
};

TEST_F(HierarchicalSystemTest, HierarchicalSystem) {
  // i, and Ω are 0 by default.
  KeplerianElements<World> elements;
  elements.eccentricity = 0;
  elements.argument_of_periapsis = 0 * Radian;
  elements.mean_anomaly = 0 * Radian;

  // Invariant: |body_indices[bodies[i]] == i| for all |i|.
  std::map<not_null<MassiveBody const*>, int> body_indices;
  std::vector<not_null<MassiveBody const*>> bodies;

  auto const new_body = [&body_indices, &bodies](Mass const& mass) {
    auto body = make_not_null_unique<MassiveBody>(mass);
    bodies.emplace_back(body.get());
    body_indices[body.get()] = body_indices.size();
    return body;
  };

  // We construct a system as follows, where the |body_indices| are, from left
  // to right, 0, 2, 3, 1.
  // |<1 m>|     |<1 m>|
  // 2     1     1     2
  //   |<   7/3 m   >|

  bodies.reserve(4);  // We'll call new_body 4 times.
  HierarchicalSystem<World> system(new_body(2 * Kilogram));
  elements.semimajor_axis = 7.0 / 3.0 * Metre;
  system.Add(new_body(2 * Kilogram), /*parent=*/bodies[0], elements);
  elements.semimajor_axis = 1 * Metre;
  system.Add(new_body(1 * Kilogram), /*parent=*/bodies[0], elements);
  elements.mean_anomaly = π * Radian;
  system.Add(new_body(1 * Kilogram), /*parent=*/bodies[1], elements);

  auto const barycentric_system = system.ConsumeBarycentricSystem();
  // primary, closest secondary, furthest secondary, child of furthest
  // secondary.
  std::vector<int> expected_order = {0, 2, 1, 3};
  for (int i = 0; i < barycentric_system.bodies.size(); ++i) {
    EXPECT_TRUE(bodies[expected_order[i]] == barycentric_system.bodies[i].get())
        << i;
  }
  std::vector<Length> x_positions;
  std::transform(barycentric_system.degrees_of_freedom.begin(),
                 barycentric_system.degrees_of_freedom.end(),
                 std::back_inserter(x_positions),
                 [](DegreesOfFreedom<World> const& dof) {
                   return (dof.position() - World::origin).coordinates().x;
                 });
  EXPECT_THAT(x_positions,
              ElementsAre(AlmostEquals(-1.5 * Metre, 0),
                          AlmostEquals(-0.5 * Metre, 0),
                          AlmostEquals(1.5 * Metre, 2),
                          AlmostEquals(0.5 * Metre, 3)));
}

TEST_F(HierarchicalSystemTest, FromMeanMotions) {
  // i, and Ω are 0 by default.
  KeplerianElements<World> elements;
  elements.eccentricity = 0;
  elements.argument_of_periapsis = 0 * Radian;
  elements.mean_anomaly = 0 * Radian;


  // Invariant: |body_indices[bodies[i]] == i| for all |i|.
  std::map<not_null<MassiveBody const*>, int> body_indices;
  std::vector<not_null<MassiveBody const*>> bodies;

  auto const new_body = [&body_indices, &bodies]() {
    auto body =
        make_not_null_unique<MassiveBody>(SIUnit<GravitationalParameter>());
    bodies.emplace_back(body.get());
    body_indices[body.get()] = body_indices.size();
    return body;
  };

  // We construct a system as follows, where the |body_indices| are, from left
  // to right, 0, 2, 1.  All bodies have unit gravitational parameter.
  // |<1 m>|
  // .     .     .
  //    |<1.5 m >|

  bodies.reserve(3);  // We'll call new_body 3 times.
  HierarchicalSystem<World> system(new_body());
  elements.mean_motion = Sqrt(3 / Pow<3>(1.5)) * Radian / Second;
  system.Add(new_body(), /*parent=*/bodies[0], elements);
  elements.mean_motion = Sqrt(2) * Radian / Second;
  system.Add(new_body(), /*parent=*/bodies[0], elements);

  auto const barycentric_system = system.ConsumeBarycentricSystem();

  std::vector<int> expected_order = {0, 2, 1};
  for (int i = 0; i < barycentric_system.bodies.size(); ++i) {
    EXPECT_TRUE(bodies[expected_order[i]] == barycentric_system.bodies[i].get())
        << i;
  }
  std::vector<Length> x_positions;
  std::transform(barycentric_system.degrees_of_freedom.begin(),
                 barycentric_system.degrees_of_freedom.end(),
                 std::back_inserter(x_positions),
                 [](DegreesOfFreedom<World> const& dof) {
                   return (dof.position() - World::origin).coordinates().x;
                 });
  EXPECT_THAT(x_positions,
              ElementsAre(AlmostEquals(-1 * Metre, 0, 1),
                          VanishesBefore(1 * Metre, 0, 1),
                          AlmostEquals(1 * Metre, 0, 1)));
}

}  // namespace internal_hierarchical_system
}  // namespace physics
}  // namespace principia
