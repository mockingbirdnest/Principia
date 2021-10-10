#include "physics/discrete_trajectory_segment.hpp"

#include <memory>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/discrete_trajectory_types.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using base::check_not_null;
using geometry::Frame;
using geometry::Instant;
using quantities::si::Second;

class DiscreteTrajectorySegmentTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  DiscreteTrajectorySegmentTest() : segments_(1) {
    auto const it = segments_.begin();
    *it = std::make_unique<DiscreteTrajectorySegment<World>>(
        DiscreteTrajectorySegmentIterator<World>(check_not_null(&segments_),
                                                 it));
    segment_ = segments_.cbegin()->get();

    segment_->Append(t0_ + 2 * Second, unmoving_origin_);
    segment_->Append(t0_ + 3 * Second, unmoving_origin_);
    segment_->Append(t0_ + 5 * Second, unmoving_origin_);
    segment_->Append(t0_ + 7 * Second, unmoving_origin_);
    segment_->Append(t0_ + 11 * Second, unmoving_origin_);
  }

  void ForgetAfter(Instant const& t) {
    segment_->ForgetAfter(t);
  }

  void ForgetBefore(Instant const& t) {
    segment_->ForgetBefore(t);
  }

  DiscreteTrajectorySegment<World>* segment_;
  internal_discrete_trajectory_types::Segments<World> segments_;
  Instant const t0_;
  DegreesOfFreedom<World> unmoving_origin_{World::origin, World::unmoving};
};

TEST_F(DiscreteTrajectorySegmentTest, Extremities) {
  {
    auto const it = segment_->begin();
    EXPECT_EQ(t0_ + 2 * Second, it->first);
  }
  {
    auto it = segment_->end();
    --it;
    EXPECT_EQ(t0_ + 11 * Second, it->first);
  }
  {
    auto const it = segment_->rbegin();
    EXPECT_EQ(t0_ + 11 * Second, it->first);
  }
  {
    auto it = segment_->rend();
    --it;
    EXPECT_EQ(t0_ + 2 * Second, it->first);
  }
}

TEST_F(DiscreteTrajectorySegmentTest, Find) {
  {
    auto const it = segment_->find(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 5 * Second, it->first);
  }
  {
    auto const it = segment_->find(t0_ + 6 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
}

TEST_F(DiscreteTrajectorySegmentTest, LowerBoundUpperBound) {
  {
    auto const it = segment_->lower_bound(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 5 * Second, it->first);
  }
  {
    auto const it = segment_->lower_bound(t0_ + 6 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->first);
  }
  {
    auto const it = segment_->lower_bound(t0_ + 12 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
  {
    auto const it = segment_->upper_bound(t0_ + 5 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->first);
  }
  {
    auto const it = segment_->upper_bound(t0_ + 6 * Second);
    EXPECT_EQ(t0_ + 7 * Second, it->first);
  }
  {
    auto const it = segment_->upper_bound(t0_ + 11 * Second);
    EXPECT_TRUE(it == segment_->end());
  }
}

TEST_F(DiscreteTrajectorySegmentTest, EmptySize) {
  EXPECT_FALSE(segment_->empty());
  EXPECT_EQ(5, segment_->size());
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterExisting) {
  ForgetAfter(t0_ + 5 * Second);
  EXPECT_EQ(t0_ + 3 * Second, segment_->rbegin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterNonexisting) {
  ForgetAfter(t0_ + 6 * Second);
  EXPECT_EQ(t0_ + 5 * Second, segment_->rbegin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetAfterPastTheEnd) {
  ForgetAfter(t0_ + 29 * Second);
  EXPECT_EQ(t0_ + 11 * Second, segment_->rbegin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeExisting) {
  ForgetBefore(t0_ + 7 * Second);
  EXPECT_EQ(t0_ + 7 * Second, segment_->begin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeNonexisting) {
  ForgetBefore(t0_ + 6 * Second);
  EXPECT_EQ(t0_ + 7 * Second, segment_->begin()->first);
}

TEST_F(DiscreteTrajectorySegmentTest, ForgetBeforeTheBeginning) {
  ForgetBefore(t0_ + 1 * Second);
  EXPECT_EQ(t0_ + 2 * Second, segment_->begin()->first);
}

}  // namespace physics
}  // namespace principia
