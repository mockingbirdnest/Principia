#include "physics/discrete_trajectory_view.hpp"

#include <ranges>

#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
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
using namespace principia::physics::_discrete_trajectory_view;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_discrete_trajectory_factories;
using namespace principia::testing_utilities::_matchers;

class DiscreteTrajectoryViewTest : public ::testing::Test {
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

TEST_F(DiscreteTrajectoryViewTest, ViewOfEntireNonemptyTrajectory) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> const view(&trajectory);

  EXPECT_EQ(trajectory.front().time, view.front().time);
  EXPECT_EQ(trajectory.front().degrees_of_freedom,
            view.front().degrees_of_freedom);
  EXPECT_EQ(trajectory.back().time, view.back().time);
  EXPECT_EQ(trajectory.back().degrees_of_freedom,
            view.back().degrees_of_freedom);

  EXPECT_EQ(trajectory.begin(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());
  EXPECT_EQ(trajectory.rbegin(), view.rbegin());
  EXPECT_EQ(trajectory.rend(), view.rend());

  EXPECT_EQ(15, view.size());

  EXPECT_EQ(t0_ + 6 * Second, view.find(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 6 * Second, view.lower_bound(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 7 * Second, view.lower_bound(t0_ + 6.1 * Second)->time);
  EXPECT_EQ(view.begin(), view.lower_bound(t0_ - 3 * Second));
  EXPECT_EQ(view.end(), view.lower_bound(t0_ + 14.1 * Second));
  EXPECT_EQ(t0_ + 7 * Second, view.upper_bound(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 7 * Second, view.upper_bound(t0_ + 6.1 * Second)->time);
  EXPECT_EQ(view.begin(), view.upper_bound(t0_ - 3 * Second));
  EXPECT_EQ(view.end(), view.upper_bound(t0_ + 14 * Second));

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 6.2 * Second),
            view.EvaluatePosition(t0_ + 6.2 * Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 6.3 * Second),
            view.EvaluateVelocity(t0_ + 6.3 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewOfEntireEmptyTrajectory) {
  DiscreteTrajectory<World> trajectory;
  DiscreteTrajectoryView<World> const view(&trajectory);

  EXPECT_EQ(trajectory.begin(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());

  EXPECT_EQ(0, view.size());

  EXPECT_EQ(view.end(), view.find(t0_ + 6 * Second));
  EXPECT_EQ(view.end(), view.lower_bound(t0_ + 6 * Second));
  EXPECT_EQ(view.end(), view.upper_bound(t0_ + 6 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByIteratorsNonempty) {
  auto const trajectory = MakeTrajectory();
  auto const it1 = std::next(trajectory.begin());
  auto const it2 = std::prev(trajectory.end());
  DiscreteTrajectoryView<World> view(&trajectory, it1, it2);

  EXPECT_EQ(it1->time, view.front().time);
  EXPECT_EQ(it1->degrees_of_freedom, view.front().degrees_of_freedom);
  EXPECT_EQ(std::prev(it2)->time, view.back().time);
  EXPECT_EQ(std::prev(it2)->degrees_of_freedom, view.back().degrees_of_freedom);

  EXPECT_EQ(it1, view.begin());
  EXPECT_EQ(it2, view.end());

  EXPECT_EQ(13, view.size());

  EXPECT_EQ(t0_ + 6 * Second, view.find(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 6 * Second, view.lower_bound(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 7 * Second, view.lower_bound(t0_ + 6.1 * Second)->time);
  EXPECT_EQ(view.begin(), view.lower_bound(t0_));
  EXPECT_EQ(view.end(), view.lower_bound(t0_ + 13.1 * Second));
  EXPECT_EQ(t0_ + 7 * Second, view.upper_bound(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 7 * Second, view.upper_bound(t0_ + 6.1 * Second)->time);
  EXPECT_EQ(view.begin(), view.upper_bound(t0_));
  EXPECT_EQ(view.end(), view.upper_bound(t0_ + 13 * Second));

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 6.2 * Second),
            view.EvaluatePosition(t0_ + 6.2 * Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 6.3 * Second),
            view.EvaluateVelocity(t0_ + 6.3 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByIteratorsEmpty) {
  auto const trajectory = MakeTrajectory();
  auto const it = std::next(trajectory.begin());
  DiscreteTrajectoryView<World> const view(&trajectory, it, it);

  EXPECT_EQ(trajectory.end(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());

  EXPECT_EQ(0, view.size());

  EXPECT_EQ(view.end(), view.find(t0_ + 6 * Second));
  EXPECT_EQ(view.end(), view.lower_bound(t0_ + 6 * Second));
  EXPECT_EQ(view.end(), view.upper_bound(t0_ + 6 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByTimesNonempty) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> const view(
      &trajectory, t0_ + 1.5 * Second, t0_ + 13.2 * Second);

  EXPECT_EQ(t0_ + 2 * Second, view.front().time);
  EXPECT_EQ(t0_ + 13 * Second, view.back().time);

  EXPECT_EQ(std::next(std::next(trajectory.begin())), view.begin());
  EXPECT_EQ(std::prev(trajectory.end()), view.end());

  EXPECT_EQ(12, view.size());

  EXPECT_EQ(t0_ + 6 * Second, view.find(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 6 * Second, view.lower_bound(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 7 * Second, view.lower_bound(t0_ + 6.1 * Second)->time);
  EXPECT_EQ(view.begin(), view.lower_bound(t0_ + 1 * Second));
  EXPECT_EQ(view.end(), view.lower_bound(t0_ + 13.1 * Second));
  EXPECT_EQ(t0_ + 7 * Second, view.upper_bound(t0_ + 6 * Second)->time);
  EXPECT_EQ(t0_ + 7 * Second, view.upper_bound(t0_ + 6.1 * Second)->time);
  EXPECT_EQ(view.begin(), view.upper_bound(t0_ + 1.2 * Second));
  EXPECT_EQ(view.end(), view.upper_bound(t0_ + 13 * Second));

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 6.2 * Second),
            view.EvaluatePosition(t0_ + 6.2 * Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 6.3 * Second),
            view.EvaluateVelocity(t0_ + 6.3 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByTimesNonemptySize0) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> const view(
      &trajectory, t0_ + 1.5 * Second, t0_ + 1.8 * Second);

  EXPECT_EQ(trajectory.end(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());

  EXPECT_FALSE(view.empty());
  EXPECT_EQ(0, view.size());

  EXPECT_EQ(view.end(), view.find(t0_ + 1.6 * Second));
  EXPECT_EQ(view.end(), view.lower_bound(t0_ + 1.6 * Second));
  EXPECT_EQ(view.end(), view.upper_bound(t0_ + 1.6 * Second));

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 1.6 * Second),
            view.EvaluatePosition(t0_ + 1.6* Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 1.7 * Second),
            view.EvaluateVelocity(t0_ + 1.7 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 1.5 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 1.5 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByTimesEmpty) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> const view(
      &trajectory, t0_ + 3.1 * Second, t0_ + 3.6 * Second);

  EXPECT_EQ(trajectory.end(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());

  EXPECT_EQ(0, view.size());

  EXPECT_EQ(view.end(), view.find(t0_ + 3.5 * Second));
  EXPECT_EQ(view.end(), view.lower_bound(t0_ + 3.5 * Second));
  EXPECT_EQ(view.end(), view.upper_bound(t0_ + 3.5 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, Deduction) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView const view1(&trajectory);
  DiscreteTrajectoryView const view2(
      &trajectory, trajectory.begin(), trajectory.end());
  DiscreteTrajectoryView const view3(&trajectory, t0_, t0_ + 14 * Second);

  EXPECT_EQ(15, view1.size());
  EXPECT_EQ(15, view2.size());
  EXPECT_EQ(15, view3.size());
}

TEST_F(DiscreteTrajectoryViewTest, Ranges) {
  auto const trajectory = MakeTrajectory();
  // Check that this compiles.
  auto it = std::ranges::begin(DiscreteTrajectoryView(&trajectory));
}

TEST_F(DiscreteTrajectoryViewTest, Restrict) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> view(&trajectory);

  EXPECT_EQ(trajectory.begin(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());
  EXPECT_EQ(t0_, trajectory.t_min());
  EXPECT_EQ(t0_ + 14 * Second, trajectory.t_max());

  EXPECT_FALSE(view.empty());
  EXPECT_EQ(15, view.size());

  view.Restrict(std::next(trajectory.begin()), trajectory.end());

  EXPECT_EQ(std::next(trajectory.begin()), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());
  EXPECT_EQ(t0_ + 1 * Second, trajectory.t_min());
  EXPECT_EQ(t0_ + 14 * Second, trajectory.t_max());

  EXPECT_FALSE(view.empty());
  EXPECT_EQ(14, view.size());

  view.Restrict(t0_ + 2.3 * Second, t0_ + 12.2 * Second);

  EXPECT_EQ(t0_ + 3 * Second, view.front().time);
  EXPECT_EQ(t0_ + 12 * Second, view.back().time);
  EXPECT_EQ(t0_ + 2.3 * Second, trajectory.t_min());
  EXPECT_EQ(t0_ + 12.2 * Second, trajectory.t_max());

  EXPECT_FALSE(view.empty());
  EXPECT_EQ(10, view.size());

  view.Restrict(t0_ + 3.2 * Second, t0_ + 3.5 * Second);

  EXPECT_EQ(trajectory.end(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());
  EXPECT_EQ(t0_ + 3.2 * Second, trajectory.t_min());
  EXPECT_EQ(t0_ + 3.5 * Second, trajectory.t_max());

  EXPECT_FALSE(view.empty());
  EXPECT_EQ(0, view.size());

  view.Restrict(trajectory.end(), trajectory.begin());

  EXPECT_EQ(trajectory.end(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());
  EXPECT_EQ(InfiniteFuture, trajectory.t_min());
  EXPECT_EQ(InfinitePast, trajectory.t_max());

  EXPECT_TRUE(view.empty());
  EXPECT_EQ(0, view.size());
}

}  // namespace physics
}  // namespace principia
