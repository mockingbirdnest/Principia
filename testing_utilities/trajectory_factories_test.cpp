#pragma once

#include "testing_utilities/trajectory_factories.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

class TrajectoryFactoriesTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
};

TEST_F(TrajectoryFactoriesTest, MakeCircularTrajectory) {
  auto const trajectory = MakeCircularTrajectory<World>(
      /*ω=*/3 * Radian / Second,
      /*r=*/2 * Metre,
      /*Δt=*/0.1 * Second, /*t1=*/
      Instant{} + 4 * Second,
      /*t2=*/Instant{} + 42 * Second);

  for (auto const& [time, degrees_of_freedom] : *trajectory) {
    Position<World> const& position = degrees_of_freedom.position();
    Velocity<World> const& velocity = degrees_of_freedom.velocity();

    EXPECT_THAT((position - World::origin).Norm(),
                AlmostEquals(2 * Metre, 0, 1));
  }
}

}  // namespace testing_utilities
}  // namespace principia
