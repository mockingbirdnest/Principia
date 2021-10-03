#include "physics/discrete_trajectory_iterator.hpp"

#include <memory>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/fake_discrete_trajectory.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using geometry::Frame;
using geometry::Instant;
using physics::DegreesOfFreedom;
using quantities::Time;
using quantities::si::Second;

class DiscreteTrajectoryIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;
  using Timeline = internal_discrete_trajectory_types::Timeline<World>;

  void FillSegment(Segments::iterator const it, Timeline const& timeline) {
    **it = DiscreteTrajectorySegment<World>{
        .self = DiscreteTrajectorySegmentIterator<World>(it),
        .timeline = timeline};
  }
  DiscreteTrajectoryIterator<World> MakeIterator(
      Segments::const_iterator const segment,
      Timeline::const_iterator const point) {
    return DiscreteTrajectoryIterator<World>(
        DiscreteTrajectorySegmentIterator<World>(segment), point);
  }


  Timeline::value_type MakeTimelineValueType(Time const& t) {
    static const DegreesOfFreedom<World> unmoving_origin(World::origin,
                                                         World::unmoving);
    return {t0_ + t, unmoving_origin};
  }

  Instant const t0_;
};

TEST_F(DiscreteTrajectoryIteratorTest, Basic) {
  DiscreteTrajectory<World> trajectory;
  for (int i = 0; i < 3; ++i) {
    trajectory.segments.push_back(
        std::make_unique<DiscreteTrajectorySegment<World>>());
  }

  auto it = trajectory.segments.begin();
  FillSegment(it++, Timeline{MakeTimelineValueType(2 * Second),
                             MakeTimelineValueType(3 * Second),
                             MakeTimelineValueType(5 * Second),
                             MakeTimelineValueType(7 * Second),
                             MakeTimelineValueType(11 * Second)});
  FillSegment(it++, Timeline{MakeTimelineValueType(13 * Second)});
  FillSegment(it++, Timeline{MakeTimelineValueType(17 * Second),
                             MakeTimelineValueType(19 * Second),
                             MakeTimelineValueType(23 * Second)});

  {
    auto segment = trajectory.segments.begin();
    auto iterator = MakeIterator(segment, (*segment)->timeline.begin());
    EXPECT_EQ(t0_ + 2 * Second, iterator->first);
    auto const current = ++iterator;
    EXPECT_EQ(t0_ + 3 * Second, iterator->first);
    EXPECT_EQ(t0_ + 3 * Second, current->first);
    auto const previous = iterator++;
    EXPECT_EQ(t0_ + 5 * Second, iterator->first);
    EXPECT_EQ(t0_ + 3 * Second, previous->first);
  }
  {
    auto segment = --trajectory.segments.end();
    auto iterator = MakeIterator(segment, (*segment)->timeline.end());
    --iterator;
    EXPECT_EQ(t0_ + 23 * Second, (*iterator).first);
    auto const current = --iterator;
    EXPECT_EQ(t0_ + 19 * Second, (*iterator).first);
    EXPECT_EQ(t0_ + 19 * Second, (*current).first);
    auto const previous = iterator--;
    EXPECT_EQ(t0_ + 17 * Second, (*iterator).first);
    EXPECT_EQ(t0_ + 19 * Second, (*previous).first);
  }
  {
    auto segment = trajectory.segments.begin();
    auto iterator = MakeIterator(segment, (*segment)->timeline.begin());
    for (int i = 0; i < 4; ++i) {
      ++iterator;
    }
    EXPECT_EQ(t0_ + 11 * Second, iterator->first);
    ++iterator;
    EXPECT_EQ(t0_ + 13 * Second, iterator->first);
    ++iterator;
    EXPECT_EQ(t0_ + 17 * Second, iterator->first);
  }
  {
    auto segment = --trajectory.segments.end();
    auto iterator = MakeIterator(segment, (*segment)->timeline.end());
    for (int i = 0; i < 3; ++i) {
      --iterator;
    }
    EXPECT_EQ(t0_ + 17 * Second, (*iterator).first);
    --iterator;
    EXPECT_EQ(t0_ + 13 * Second, (*iterator).first);
    --iterator;
    EXPECT_EQ(t0_ + 11 * Second, (*iterator).first);
  }
}

}  // namespace physics
}  // namespace principia
