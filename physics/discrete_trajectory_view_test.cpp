#include "physics/discrete_trajectory_view.hpp"

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

TEST_F(DiscreteTrajectoryViewTest, ViewOfEntireTrajectory) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> view(&trajectory);

  EXPECT_EQ(trajectory.front().time, view.front().time);
  EXPECT_EQ(trajectory.front().degrees_of_freedom,
            view.front().degrees_of_freedom);
  EXPECT_EQ(trajectory.back().time, view.back().time);
  EXPECT_EQ(trajectory.back().degrees_of_freedom,
            view.back().degrees_of_freedom);

  EXPECT_EQ(trajectory.begin(), view.begin());
  EXPECT_EQ(trajectory.end(), view.end());

  EXPECT_FALSE(view.empty());
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

  EXPECT_EQ(trajectory.t_min(), view.t_min());
  EXPECT_EQ(trajectory.t_max(), view.t_max());

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 6.2 * Second),
            view.EvaluatePosition(t0_ + 6.2 * Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 6.3 * Second),
            view.EvaluateVelocity(t0_ + 6.3 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByIteratorsNonempty) {
  auto const trajectory = MakeTrajectory();
  auto const it1 = std::next(trajectory.begin());
  auto const it2 = std::prev(trajectory.end());
  DiscreteTrajectoryView<World> view(&trajectory, it1, it2);

  EXPECT_EQ(it1->time, view.front().time);
  EXPECT_EQ(it1->degrees_of_freedom, view.front().degrees_of_freedom);
  EXPECT_EQ(it2->time, view.back().time);
  EXPECT_EQ(it2->degrees_of_freedom, view.back().degrees_of_freedom);

  EXPECT_EQ(it1, view.begin());
  EXPECT_EQ(it2, view.end());

  EXPECT_FALSE(view.empty());
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

  EXPECT_EQ(t0_ + 1 * Second, view.t_min());
  EXPECT_EQ(t0_ + 13 * Second, view.t_max());

  EXPECT_EQ(trajectory.EvaluatePosition(t0_ + 6.2 * Second),
            view.EvaluatePosition(t0_ + 6.2 * Second));
  EXPECT_EQ(trajectory.EvaluateVelocity(t0_ + 6.3 * Second),
            view.EvaluateVelocity(t0_ + 6.3 * Second));
  EXPECT_EQ(trajectory.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second),
            view.EvaluateDegreesOfFreedom(t0_ + 6.4 * Second));
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByIteratorsEmpty) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> view(
      &trajectory, std::prev(trajectory.end()), std::next(trajectory.begin()));
  EXPECT_EQ(InfiniteFuture, view.t_min());
  EXPECT_EQ(InfinitePast, view.t_max());
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByTimesNonempty) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> view(
      &trajectory, t0_ + 1.5 * Second, t0_ + 13.2 * Second);
  EXPECT_EQ(t0_ + 2 * Second, view.t_min());
  EXPECT_EQ(t0_ + 13 * Second, view.t_max());
}

TEST_F(DiscreteTrajectoryViewTest, ViewDefinedByTimesEmpty) {
  auto const trajectory = MakeTrajectory();
  DiscreteTrajectoryView<World> view(
      &trajectory, t0_ + 13 * Second, t0_ + 3 * Second);
  EXPECT_EQ(InfiniteFuture, view.t_min());
  EXPECT_EQ(InfinitePast, view.t_max());
}

}  // namespace physics
}  // namespace principia
