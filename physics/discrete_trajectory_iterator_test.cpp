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
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using base::check_not_null;
using base::make_not_null_unique;
using base::not_null;
using geometry::Frame;
using geometry::Instant;
using physics::DegreesOfFreedom;
using quantities::Time;
using quantities::si::Second;
using ::testing::Return;

class DiscreteTrajectoryIteratorTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
  using Segments = internal_discrete_trajectory_types::Segments<World>;

  DiscreteTrajectoryIteratorTest()
      : segments_(MakeSegments(3)) {
    auto it = segments_->begin();
    {
      auto& segment1 = *it;
      segment1.Append(t0_ + 2 * Second, unmoving_origin_);
      segment1.Append(t0_ + 3 * Second, unmoving_origin_);
      segment1.Append(t0_ + 5 * Second, unmoving_origin_);
      segment1.Append(t0_ + 7 * Second, unmoving_origin_);
      segment1.Append(t0_ + 11 * Second, unmoving_origin_);
    }

    ++it;
    {
      auto& segment2 = *it;
      segment2.Append(t0_ + 13 * Second, unmoving_origin_);
    }

    ++it;
    {
      auto& segment3 = *it;
      segment3.Append(t0_ + 13 * Second, unmoving_origin_);
      segment3.Append(t0_ + 17 * Second, unmoving_origin_);
      segment3.Append(t0_ + 19 * Second, unmoving_origin_);
      segment3.Append(t0_ + 23 * Second, unmoving_origin_);
    }
  }

  DiscreteTrajectoryIterator<World> MakeBegin(
      Segments::const_iterator const it) {
    return DiscreteTrajectoryIterator<World>(
        DiscreteTrajectorySegmentIterator<World>(segments_.get(), it),
        it->timeline_begin());
  }

  DiscreteTrajectoryIterator<World> MakeEnd(Segments::const_iterator it) {
    return DiscreteTrajectoryIterator<World>(
        DiscreteTrajectorySegmentIterator<World>(segments_.get(), ++it),
        std::nullopt);
  }

  // Constructs a list of |n| segments which are properly initialized.
  // TODO(phl): Move to a central place.
  static not_null<std::unique_ptr<Segments>> MakeSegments(const int n) {
    auto segments = make_not_null_unique<Segments>(n);
    for (auto it = segments->begin(); it != segments->end(); ++it) {
      *it = DiscreteTrajectorySegment<World>(
          DiscreteTrajectorySegmentIterator<World>(segments.get(), it));
    }
    return segments;
  }

  not_null<std::unique_ptr<Segments>> segments_;
  Instant const t0_;
  DegreesOfFreedom<World> const unmoving_origin_{World::origin,
                                                 World::unmoving};
};

TEST_F(DiscreteTrajectoryIteratorTest, ForwardOneSegment) {
  auto segment = segments_->begin();
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
  auto segment = --segments_->end();
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
  auto segment = segments_->begin();
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
  auto segment = --segments_->end();
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
  auto it1 = MakeEnd(--segments_->end());
  for (int i = 0; i < 3; ++i) {
    --it1;
  }
  EXPECT_EQ(t0_ + 17 * Second, (*it1).first);
  --it1;
  EXPECT_EQ(t0_ + 13 * Second, (*it1).first);

  auto it2 = MakeBegin(segments_->begin());
  for (int i = 0; i < 4; ++i) {
    ++it2;
  }
  EXPECT_EQ(t0_ + 11 * Second, it2->first);
  ++it2;
  EXPECT_EQ(t0_ + 13 * Second, it2->first);

  EXPECT_EQ(it1, it2);
  EXPECT_NE(it1, MakeBegin(segments_->begin()));
  EXPECT_NE(it2, MakeEnd(--segments_->end()));
  EXPECT_NE(MakeBegin(segments_->begin()), MakeEnd(--segments_->end()));
}

}  // namespace physics
}  // namespace principia
