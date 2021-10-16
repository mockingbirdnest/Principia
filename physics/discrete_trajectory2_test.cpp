#include "physics/discrete_trajectory2.hpp"

#include <vector>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"

namespace principia {
namespace physics {

using geometry::Displacement;
using geometry::Frame;
using geometry::Instant;
using geometry::Velocity;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::Componentwise;
using testing_utilities::NewLinearTrajectoryTimeline;
using ::testing::ElementsAre;
using ::testing::Eq;

class DiscreteTrajectory2Test : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  // Constructs a trajectory with three 5-second segments.
  DiscreteTrajectory2<World> MakeTrajectory() {
    DiscreteTrajectory2<World> trajectory;
    std::optional<DegreesOfFreedom<World>> last_degrees_of_freedom;

    Velocity<World> const v1({1 * Metre / Second,
                              0 * Metre / Second,
                              0 * Metre / Second});
    for (auto const& [t, degrees_of_freedom] :
         NewLinearTrajectoryTimeline(v1,
                                     /*Δt=*/1 * Second,
                                     /*t1=*/t0_,
                                     /*t2=*/t0_ + 5 * Second)) {
      last_degrees_of_freedom = degrees_of_freedom;
      trajectory.Append(t, degrees_of_freedom);
    }

    trajectory.NewSegment();
    Velocity<World> const v2({0 * Metre / Second,
                              1 * Metre / Second,
                              0 * Metre / Second});
    for (auto const& [t, degrees_of_freedom] :
        NewLinearTrajectoryTimeline(DegreesOfFreedom<World>(
                                        last_degrees_of_freedom->position(),
                                        v2),
                                     /*Δt=*/1 * Second,
                                     /*t1=*/t0_ + 5 * Second,
                                     /*t2=*/t0_ + 10 * Second)) {
      last_degrees_of_freedom = degrees_of_freedom;
      trajectory.Append(t, degrees_of_freedom);
    }

    trajectory.NewSegment();
    Velocity<World> const v3({0 * Metre / Second,
                              0 * Metre / Second,
                              1 * Metre / Second});
    for (auto const& [t, degrees_of_freedom] :
        NewLinearTrajectoryTimeline(DegreesOfFreedom<World>(
                                        last_degrees_of_freedom->position(),
                                        v3),
                                     /*Δt=*/1 * Second,
                                     /*t1=*/t0_ + 10 * Second,
                                     /*t2=*/t0_ + 15 * Second)) {
      trajectory.Append(t, degrees_of_freedom);
    }

    return trajectory;
  }

  Instant const t0_;
};

TEST_F(DiscreteTrajectory2Test, Make) {
  auto const trajectory = MakeTrajectory();
}

TEST_F(DiscreteTrajectory2Test, IterateForward) {
  auto const trajectory = MakeTrajectory();
  std::vector<Instant> times;
  for (auto const& [t, _] : trajectory) {
    times.push_back(t);
  }
  EXPECT_THAT(times,
              ElementsAre(t0_,
                          t0_ + 1 * Second,
                          t0_+ 2 * Second,
                          t0_+ 3 * Second,
                          t0_+ 4 * Second,
                          t0_+ 5 * Second,
                          t0_+ 6 * Second,
                          t0_+ 7 * Second,
                          t0_+ 8 * Second,
                          t0_+ 9 * Second,
                          t0_+ 10 * Second,
                          t0_+ 11 * Second,
                          t0_+ 12 * Second,
                          t0_+ 13 * Second,
                          t0_+ 14 * Second));
}

TEST_F(DiscreteTrajectory2Test, IterateBackward) {
  auto const trajectory = MakeTrajectory();
  std::vector<Instant> times;
  for (auto it = trajectory.rbegin(); it != trajectory.rend(); ++it) {
    times.push_back(it->first);
  }
  EXPECT_THAT(times,
              ElementsAre(t0_ + 14 * Second,
                          t0_ + 13 * Second,
                          t0_ + 12 * Second,
                          t0_ + 11 * Second,
                          t0_ + 10 * Second,
                          t0_ + 9 * Second,
                          t0_ + 8 * Second,
                          t0_ + 7 * Second,
                          t0_ + 6 * Second,
                          t0_ + 5 * Second,
                          t0_ + 4 * Second,
                          t0_ + 3 * Second,
                          t0_ + 2 * Second,
                          t0_ + 1 * Second,
                          t0_));
}

TEST_F(DiscreteTrajectory2Test, Empty) {
  DiscreteTrajectory2<World> trajectory;
  EXPECT_TRUE(trajectory.empty());
  trajectory = MakeTrajectory();
  EXPECT_FALSE(trajectory.empty());
}

TEST_F(DiscreteTrajectory2Test, Size) {
  DiscreteTrajectory2<World> trajectory;
  EXPECT_EQ(0, trajectory.size());
  trajectory = MakeTrajectory();
  EXPECT_EQ(15, trajectory.size());
}

TEST_F(DiscreteTrajectory2Test, Find) {
  auto const trajectory = MakeTrajectory();
  {
    auto const it = trajectory.find(t0_ + 3 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 3 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({3 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}));
  }
  {
    auto const it = trajectory.find(t0_ + 13 * Second);
    auto const& [t, degrees_of_freedom] = *it;
    EXPECT_EQ(t, t0_ + 13 * Second);
    EXPECT_EQ(degrees_of_freedom.position(),
              World::origin + Displacement<World>({4 * Metre,
                                                   4 * Metre,
                                                   3 * Metre}));
  }
  {
    auto const it = trajectory.find(t0_ + 3.14 * Second);
    EXPECT_TRUE(it == trajectory.end());
  }
}

}  // namespace physics
}  // namespace principia
