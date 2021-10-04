#include "physics/discrete_trajectory_iterator.hpp"

#include <memory>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/mock_discrete_trajectory_segment.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using geometry::Frame;
using geometry::Instant;
using physics::DegreesOfFreedom;
using quantities::Time;
using quantities::si::Second;
using ::testing::Return;

namespace {

//TODO(phl):comment.
template<typename Frame>
class FakeDiscreteTrajectorySegment
    : public MockDiscreteTrajectorySegment<Frame> {
 public:
  FakeDiscreteTrajectorySegment() = default;

  DiscreteTrajectorySegmentIterator<Frame> self;
  internal_discrete_trajectory_types::Timeline<Frame> timeline;
};

}  // namespace

class DiscreteTrajectoryIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;
  using Timeline = internal_discrete_trajectory_types::Timeline<World>;

  FakeDiscreteTrajectorySegment<World>* DownCast(
      std::unique_ptr<DiscreteTrajectorySegment<World>> const& segment) {
    return dynamic_cast<FakeDiscreteTrajectorySegment<World>*>(segment.get());
  }

  void FillSegment(Segments::iterator const it, Timeline const& timeline) {
    auto* const segment = DownCast(*it);
    segment->self = DiscreteTrajectorySegmentIterator<World>(it);
    segment->timeline = timeline;
  }

  DiscreteTrajectoryIterator<World> MakeBegin(
      Segments::const_iterator const it) {
    return DiscreteTrajectoryIterator<World>(
        DiscreteTrajectorySegmentIterator<World>(it),
        DownCast(*it)->timeline.begin());
  }

  DiscreteTrajectoryIterator<World> MakeEnd(
      Segments::const_iterator const it) {
    return DiscreteTrajectoryIterator<World>(
        DiscreteTrajectorySegmentIterator<World>(it),
        DownCast(*it)->timeline.end());
  }

  Timeline::value_type MakeTimelineValueType(Time const& t) {
    static const DegreesOfFreedom<World> unmoving_origin(World::origin,
                                                         World::unmoving);
    return {t0_ + t, unmoving_origin};
  }

  Instant const t0_;
};

TEST_F(DiscreteTrajectoryIteratorTest, Basic) {
  auto owned_mock1 = std::make_unique<FakeDiscreteTrajectorySegment<World>>();
  auto owned_mock2 = std::make_unique<FakeDiscreteTrajectorySegment<World>>();
  auto owned_mock3 = std::make_unique<FakeDiscreteTrajectorySegment<World>>();
  auto const& mock1 = *owned_mock1;
  auto const& mock2 = *owned_mock2;
  auto const& mock3 = *owned_mock3;

  Segments segments;
  segments.push_back(std::move(owned_mock1));
  auto const it1 = --segments.end();
  segments.push_back(std::move(owned_mock2));
  auto const it2 = --segments.end();
  segments.push_back(std::move(owned_mock3));
  auto const it3 = --segments.end();

  FillSegment(it1, Timeline{MakeTimelineValueType(2 * Second),
                            MakeTimelineValueType(3 * Second),
                            MakeTimelineValueType(5 * Second),
                            MakeTimelineValueType(7 * Second),
                            MakeTimelineValueType(11 * Second)});
  FillSegment(it2, Timeline{MakeTimelineValueType(13 * Second)});
  FillSegment(it3, Timeline{MakeTimelineValueType(13 * Second),  // Duplicated.
                            MakeTimelineValueType(17 * Second),
                            MakeTimelineValueType(19 * Second),
                            MakeTimelineValueType(23 * Second)});

  // This must happen *after* the segments have been set up.
  EXPECT_CALL(mock1, begin()).WillRepeatedly(Return(MakeBegin(it1)));
  EXPECT_CALL(mock1, end()).WillRepeatedly(Return(MakeEnd(it1)));
  EXPECT_CALL(mock2, begin()).WillRepeatedly(Return(MakeBegin(it2)));
  EXPECT_CALL(mock2, end()).WillRepeatedly(Return(MakeEnd(it2)));
  EXPECT_CALL(mock3, begin()).WillRepeatedly(Return(MakeBegin(it3)));
  EXPECT_CALL(mock3, end()).WillRepeatedly(Return(MakeEnd(it3)));

  // Iteration in one segment.
  {
    auto segment = segments.begin();
    auto iterator = MakeBegin(segment);
    EXPECT_EQ(t0_ + 2 * Second, iterator->first);
    auto const current = ++iterator;
    EXPECT_EQ(t0_ + 3 * Second, iterator->first);
    EXPECT_EQ(t0_ + 3 * Second, current->first);
    auto const previous = iterator++;
    EXPECT_EQ(t0_ + 5 * Second, iterator->first);
    EXPECT_EQ(t0_ + 3 * Second, previous->first);
  }
  {
    auto segment = --segments.end();
    auto iterator = MakeEnd(segment);
    --iterator;
    EXPECT_EQ(t0_ + 23 * Second, (*iterator).first);
    auto const current = --iterator;
    EXPECT_EQ(t0_ + 19 * Second, (*iterator).first);
    EXPECT_EQ(t0_ + 19 * Second, (*current).first);
    auto const previous = iterator--;
    EXPECT_EQ(t0_ + 17 * Second, (*iterator).first);
    EXPECT_EQ(t0_ + 19 * Second, (*previous).first);
  }

  // Iteration across segments.
  {
    auto segment = segments.begin();
    auto iterator = MakeBegin(segment);
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
    auto segment = --segments.end();
    auto iterator = MakeEnd(segment);
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
