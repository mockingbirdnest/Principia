#pragma once

#include "testing_utilities/discrete_trajectory_factories.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace testing_utilities {

using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::InnerProduct;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using quantities::Acceleration;
using quantities::Pow;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

class DiscreteTrajectoryFactoriesTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
};

TEST_F(DiscreteTrajectoryFactoriesTest, NewLinearTrajectoryTimeline) {
  auto const timeline = NewLinearTrajectoryTimeline<World>(
      /*degrees_of_freedom=*/
      DegreesOfFreedom<World>(
          World::origin +
              Displacement<World>({30 * Metre, 40 * Metre, 50 * Metre}),
          Velocity<World>(
              {6 * Metre / Second, 5 * Metre / Second, 4 * Metre / Second})),
      /*Δt=*/0.1 * Second,
      /*t0=*/Instant(),
      /*t1=*/Instant() + 4 * Second,
      /*t2=*/Instant() + 42 * Second);

  for (auto const& [time, degrees_of_freedom] : timeline) {
    Position<World> const& position = degrees_of_freedom.position();
    Velocity<World> const& velocity = degrees_of_freedom.velocity();

    EXPECT_THAT((position - World::origin).Norm(),
                AlmostEquals(Sqrt(5000 + 1160 * (time - Instant()) / Second +
                                  77 * Pow<2>((time - Instant()) / Second)) *
                                 Metre,
                             0, 1));
    EXPECT_THAT(velocity.Norm(), AlmostEquals(Sqrt(77) * Metre / Second, 0, 0));
  }
  EXPECT_THAT(timeline.begin()->time,
              AlmostEquals(Instant() + 4 * Second, 0));
  EXPECT_THAT(timeline.rbegin()->time,
              AlmostEquals(Instant() + 41.9 * Second, 46));
  EXPECT_EQ(380, timeline.size());
}

TEST_F(DiscreteTrajectoryFactoriesTest, NewCircularTrajectoryTimeline) {
  auto const timeline = NewCircularTrajectoryTimeline<World>(
      /*ω=*/3 * Radian / Second,
      /*r=*/2 * Metre,
      /*Δt=*/0.1 * Second,
      /*t1=*/Instant() + 4 * Second,
      /*t2=*/Instant() + 42 * Second);

  for (auto const& [time, degrees_of_freedom] : timeline) {
    Position<World> const& position = degrees_of_freedom.position();
    Velocity<World> const& velocity = degrees_of_freedom.velocity();

    EXPECT_THAT((position - World::origin).Norm(),
                AlmostEquals(2 * Metre, 0, 1));
    EXPECT_THAT(InnerProduct(position - World::origin, velocity),
                VanishesBefore(1 * Metre * Metre / Second, 0, 8));
  }
  EXPECT_THAT(timeline.begin()->time,
              AlmostEquals(Instant() + 4 * Second, 0));
  EXPECT_THAT(timeline.rbegin()->time,
              AlmostEquals(Instant() + 41.9 * Second, 46));
  EXPECT_EQ(380, timeline.size());
}

}  // namespace testing_utilities
}  // namespace principia
