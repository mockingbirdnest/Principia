#include "physics/trajectory_view.hpp"

#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/trajectory.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace physics {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_trajectory;
using namespace principia::physics::_trajectory_view;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_discrete_trajectory_factories;
using namespace principia::testing_utilities::_matchers;

class TrajectoryViewTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  // Constructs a trajectory with three 5-second segments starting at `t0` and
  // the given initial `degrees_of_freedom`.
  // TODO(phl): This code is duplicated in `DiscreteTrajectoryTest`.
  static DiscreteTrajectory<World> MakeTrajectory(
      Instant const& t0,
      DegreesOfFreedom<World> const& degrees_of_freedom) {
    DiscreteTrajectory<World> trajectory;
    std::optional<DegreesOfFreedom<World>> last_degrees_of_freedom;

    for (auto const& [t, degrees_of_freedom] :
         NewLinearTrajectoryTimeline(degrees_of_freedom,
                                     /*Δt=*/1 * Second,
                                     /*t1=*/t0,
                                     /*t2=*/t0 + 5 * Second)) {
      last_degrees_of_freedom = degrees_of_freedom;
      EXPECT_OK(trajectory.Append(t, degrees_of_freedom));
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
                                     /*t1=*/t0 + 5 * Second,
                                     /*t2=*/t0 + 10 * Second)) {
      last_degrees_of_freedom = degrees_of_freedom;
      EXPECT_OK(trajectory.Append(t, degrees_of_freedom));
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
                                     /*t1=*/t0 + 10 * Second,
                                     /*t2=*/t0 + 15 * Second)) {
      EXPECT_OK(trajectory.Append(t, degrees_of_freedom));
    }

    return trajectory;
  }

  DiscreteTrajectory<World> MakeTrajectory() {
    Velocity<World> const v1({1 * Metre / Second,
                              0 * Metre / Second,
                              0 * Metre / Second});
    return MakeTrajectory(t0_, DegreesOfFreedom<World>(World::origin, v1));
  }

  Instant const t0_;
};

TEST_F(TrajectoryViewTest, ViewOfEntireNonemptyTrajectory) {
  auto const trajectory = MakeTrajectory();
  TrajectoryView<World> view(&trajectory);

  EXPECT_FALSE(view.empty());

  EXPECT_EQ(trajectory.t_min(), view.t_min());
  EXPECT_EQ(trajectory.t_max(), view.t_max());

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 6.2 * Second),
            view.EvaluatePosition(t0_ + 6.2 * Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 6.3 * Second),
            view.EvaluateVelocity(t0_ + 6.3 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second));
}

TEST_F(TrajectoryViewTest, ViewOfEntireEmptyTrajectory) {
  DiscreteTrajectory<World> trajectory;
  TrajectoryView<World> view(&trajectory);

  EXPECT_TRUE(view.empty());

  EXPECT_EQ(InfiniteFuture, view.t_min());
  EXPECT_EQ(InfinitePast, view.t_max());
}

TEST_F(TrajectoryViewTest, ViewDefinedByTimesNonempty) {
  auto const trajectory = MakeTrajectory();
  TrajectoryView<World> view(
      &trajectory, t0_ + 1.5 * Second, t0_ + 13.2 * Second);

  EXPECT_FALSE(view.empty());

  EXPECT_EQ(t0_ + 1.5 * Second, view.t_min());
  EXPECT_EQ(t0_ + 13.2 * Second, view.t_max());

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 6.2 * Second),
            view.EvaluatePosition(t0_ + 6.2 * Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 6.3 * Second),
            view.EvaluateVelocity(t0_ + 6.3 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second));
}

TEST_F(TrajectoryViewTest, ViewDefinedByTimesEmpty) {
  auto const trajectory = MakeTrajectory();
  TrajectoryView<World> view(
      &trajectory, t0_ + 3.2 * Second, t0_ + 3.1 * Second);

  EXPECT_TRUE(view.empty());

  EXPECT_EQ(InfiniteFuture, view.t_min());
  EXPECT_EQ(InfinitePast, view.t_max());
}

}  // namespace physics
}  // namespace principia
