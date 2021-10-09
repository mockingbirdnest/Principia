#include "physics/discrete_trajectory_iterator.hpp"

#include <memory>
#include <optional>

#include "base/not_null.hpp"
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

using base::check_not_null;
using geometry::Frame;
using geometry::Instant;
using physics::DegreesOfFreedom;
using quantities::Time;
using quantities::si::Second;
using ::testing::Return;

namespace {

// A trajectory that holds a real timeline, where we mock multiple methods to
// return data from that timeline.
template<typename Frame>
class FakeDiscreteTrajectorySegment
    : public MockDiscreteTrajectorySegment<Frame> {
 public:
  FakeDiscreteTrajectorySegment() = default;

  internal_discrete_trajectory_types::Timeline<Frame> timeline;
};

}  // namespace

class DiscreteTrajectoryIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using Segments = internal_discrete_trajectory_types::Segments<World>;
  using Timeline = internal_discrete_trajectory_types::Timeline<World>;

  DiscreteTrajectoryIteratorTest() {
    // Set up a fake trajectory with 3 segments.  After construction, the mocks
    // are owned by |segments_|.
    auto owned_mock1 = std::make_unique<FakeDiscreteTrajectorySegment<World>>();
    auto owned_mock2 = std::make_unique<FakeDiscreteTrajectorySegment<World>>();
    auto owned_mock3 = std::make_unique<FakeDiscreteTrajectorySegment<World>>();
    auto const& mock1 = *owned_mock1;
    auto const& mock2 = *owned_mock2;
    auto const& mock3 = *owned_mock3;

    segments_.push_back(std::move(owned_mock1));
    auto const it1 = --segments_.end();
    segments_.push_back(std::move(owned_mock2));
    auto const it2 = --segments_.end();
    segments_.push_back(std::move(owned_mock3));
    auto const it3 = --segments_.end();

    FillSegment(it1,
                Timeline{MakeTimelineValueType(2 * Second),
                         MakeTimelineValueType(3 * Second),
                         MakeTimelineValueType(5 * Second),
                         MakeTimelineValueType(7 * Second),
                         MakeTimelineValueType(11 * Second)});
    FillSegment(it2, Timeline{MakeTimelineValueType(13 * Second)});
    FillSegment(it3,
                Timeline{MakeTimelineValueType(13 * Second),  // Duplicated.
                         MakeTimelineValueType(17 * Second),
                         MakeTimelineValueType(19 * Second),
                         MakeTimelineValueType(23 * Second)});

    // This must happen *after* the segments have been set up.
    EXPECT_CALL(mock1, timeline_begin())
        .WillRepeatedly(Return(timeline_begin(it1)));
    EXPECT_CALL(mock1, timeline_end())
        .WillRepeatedly(Return(timeline_end(it1)));
    EXPECT_CALL(mock2, timeline_begin())
        .WillRepeatedly(Return(timeline_begin(it2)));
    EXPECT_CALL(mock2, timeline_end())
        .WillRepeatedly(Return(timeline_end(it2)));
    EXPECT_CALL(mock3, timeline_begin())
        .WillRepeatedly(Return(timeline_begin(it3)));
    EXPECT_CALL(mock3, timeline_end())
        .WillRepeatedly(Return(timeline_end(it3)));
  }

  FakeDiscreteTrajectorySegment<World>* DownCast(
      std::unique_ptr<DiscreteTrajectorySegment<World>> const& segment) {
    return dynamic_cast<FakeDiscreteTrajectorySegment<World>*>(segment.get());
  }

  void FillSegment(Segments::iterator const it, Timeline const& timeline) {
    auto* const segment = DownCast(*it);
    segment->timeline = timeline;
  }

  DiscreteTrajectoryIterator<World> MakeBegin(
      Segments::const_iterator const it) {
    return DiscreteTrajectoryIterator<World>(
        DiscreteTrajectorySegmentIterator<World>(
        check_not_null(&segments_), it),
        timeline_begin(it));
  }

  DiscreteTrajectoryIterator<World> MakeEnd(Segments::const_iterator it) {
    return DiscreteTrajectoryIterator<World>(
        DiscreteTrajectorySegmentIterator<World>(check_not_null(&segments_),
                                                 ++it),
        std::nullopt);
  }

  Timeline::value_type MakeTimelineValueType(Time const& t) {
    static const DegreesOfFreedom<World> unmoving_origin(World::origin,
                                                         World::unmoving);
    return {t0_ + t, unmoving_origin};
  }

  internal_discrete_trajectory_types::Timeline<World>::const_iterator
  timeline_begin(Segments::const_iterator const it) {
    return DownCast(*it)->timeline.begin();
  }

  internal_discrete_trajectory_types::Timeline<World>::const_iterator
  timeline_end(Segments::const_iterator const it) {
    return DownCast(*it)->timeline.end();
  }

  Segments segments_;
  Instant const t0_;
};

TEST_F(DiscreteTrajectoryIteratorTest, ForwardOneSegment) {
  auto segment = segments_.begin();
  auto iterator = MakeBegin(segment);
  EXPECT_EQ(t0_ + 2 * Second, iterator->first);
  auto const current = ++iterator;
  EXPECT_EQ(t0_ + 3 * Second, iterator->first);
  EXPECT_EQ(t0_ + 3 * Second, current->first);
  auto const previous = iterator++;
  EXPECT_EQ(t0_ + 5 * Second, iterator->first);
  EXPECT_EQ(t0_ + 3 * Second, previous->first);
}

TEST_F(DiscreteTrajectoryIteratorTest, BackwardOneSegment) {
  auto segment = --segments_.end();
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

TEST_F(DiscreteTrajectoryIteratorTest, ForwardAcrossSegments) {
  auto segment = segments_.begin();
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

TEST_F(DiscreteTrajectoryIteratorTest, BackwardAcrossSegments) {
  auto segment = --segments_.end();
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

TEST_F(DiscreteTrajectoryIteratorTest, Equality) {
  // Construct two iterators that denote the time 13 * Second but in different
  // segments.
  auto it1 = MakeEnd(--segments_.end());
  for (int i = 0; i < 3; ++i) {
    --it1;
  }
  EXPECT_EQ(t0_ + 17 * Second, (*it1).first);
  --it1;
  EXPECT_EQ(t0_ + 13 * Second, (*it1).first);

  auto it2 = MakeBegin(segments_.begin());
  for (int i = 0; i < 4; ++i) {
    ++it2;
  }
  EXPECT_EQ(t0_ + 11 * Second, it2->first);
  ++it2;
  EXPECT_EQ(t0_ + 13 * Second, it2->first);

  EXPECT_EQ(it1, it2);
  EXPECT_NE(it1, MakeBegin(segments_.begin()));
  EXPECT_NE(it2, MakeEnd(--segments_.end()));
  EXPECT_NE(MakeBegin(segments_.begin()), MakeEnd(--segments_.end()));
}

}  // namespace physics
}  // namespace principia
